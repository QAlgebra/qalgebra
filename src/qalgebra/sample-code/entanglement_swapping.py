import qalgebra as qa
from sympy import sqrt
import sympy as sp
import numpy as np
from qalgebra import convert_to_sympy_matrix, KetBra
import time
hs = [qa.LocalSpace(i, dimension=2) for i in range(5)]
def hadamard_transform(hs):
    return (
        qa.LocalSigma(0, 0, hs=hs)  # LocalSigma(i,j) = |i⟩⟨j|
        + qa.LocalSigma(0, 1, hs=hs)
        + qa.LocalSigma(1, 0, hs=hs)
        - qa.LocalSigma(1, 1, hs=hs)
    ) / sqrt(2)
def ket(*indices, hs):
    local_spaces = hs.local_factors # for multiple indices, hs should be product space
    local_kets = [ls.basis_state(i) for (ls, i) in zip(local_spaces, indices)]
    return qa.TensorKet.create(*local_kets)
def ketbra(indices_ket, indices_bra, hs):
    return qa.KetBra(ket(*indices_ket, hs=hs), ket(*indices_bra, hs=hs))

# Controlled not
def cnot(hs1, hs2):
    hs = hs1*hs2
    return (
        ketbra((0,0), (0,0), hs)
        + ketbra((0,1), (0,1), hs)
        + ketbra((1,0), (1,1), hs)
        + ketbra((1, 1), (1,0), hs)
    )

# Pauli operators to reconstruct the repeated qubit and to help with entanglement
def sigmax(hs):
    return qa.LocalSigma(0, 1, hs=hs) + qa.LocalSigma(1, 0, hs=hs)
def sigmay(hs):
    return -sp.I*qa.LocalSigma(0,1,hs=hs) + sp.I*qa.LocalSigma(1,0,hs=hs)
def sigmaz(hs):
    return qa.LocalSigma(0,0, hs=hs)-qa.LocalSigma(1,1, hs=hs)

# The classic Bell states
def bell_state(state, hs):
    bs = {
        '00': ket(0, 0, hs=hs) + ket(1, 1, hs=hs),
        '01': ket(0, 0, hs=hs) - ket(1, 1, hs=hs),
        '10': ket(0, 1, hs=hs) + ket(1, 0, hs=hs),
        '11': ket(0, 1, hs=hs) - ket(1, 0, hs=hs),
    }
    return bs[state] / sqrt(2)

# Make positive operator-valued measurement projections
def proj(state):
    return KetBra(state, state).expand()

def get_state_from_pure_dm(dm, factor=None):
    dm_sympy = convert_to_sympy_matrix(dm.doit())
    triples = dm_sympy.eigenvects()
    for i, triple in enumerate(triples):
        eval1, mult, evecs = triple
        if eval1 != 0:
            assert len(evecs)==1
            x = evecs[0][:]
            if factor is not None:
                x = (factor*evecs[0])[:]
            return x

show_intermediate = False
start = time.time()
# Put the data qubit in a general superposition
a, b = sp.symbols('a,b')
initial_state = a*hs[0].basis_state(0)+b*hs[0].basis_state(1)

# The repeated data qubit
expected_final = a*hs[4].basis_state(0)+b*hs[4].basis_state(1)
dm_expected_final = KetBra(expected_final, expected_final).expand()

# The combination of the data qubit, the link qubits and the repeated qubit
q = initial_state * hs[1].basis_state(0) * hs[2].basis_state(0) *\
    hs[3].basis_state(0)* hs[4].basis_state(0)

# Entangle qubits 1 and 2 in a Bell Psi minus state
x1 = sigmax(hs[1])
x2 = sigmax(hs[2])
h = hadamard_transform(hs[1])
not1 = cnot(hs[1], hs[2])
psi_23 = (h * x2 * x1 * q)
psi_23 = (not1 * psi_23).expand()

# Entangle qubits 3 and 4 likewise
x1 = sigmax(hs[3])
x2 = sigmax(hs[4])
h = hadamard_transform(hs[3])
not1 = cnot(hs[3], hs[4])
psi_23 = (h * x2 * x1 * psi_23).expand()
psi_23 = (not1 * psi_23).expand()

# Prepare to measure qubits 2 and 3 in the Bell basis
phi_plus = bell_state('00', hs[2]*hs[3])
phi_minus = bell_state('01', hs[2]*hs[3])
psi_plus= bell_state('10', hs[2]*hs[3])
psi_minus = bell_state('11', hs[2]*hs[3])

# Prepare to measure qubits 0 and 1 in the Bell basis
phi_plus2 = bell_state('00', hs[0]*hs[1])
phi_minus2 = bell_state('01', hs[0]*hs[1])
psi_plus2= bell_state('10', hs[0]*hs[1])
psi_minus2 = bell_state('11', hs[0]*hs[1])

# Cycle through all 16 possible measurement outcomes
for i1, psi1 in enumerate([phi_plus, phi_minus, psi_plus, psi_minus]):
    
    # Measure
    measure1 = KetBra(psi1, psi1)
    psi_123 = measure1 * psi_23
    for i2, psi2 in enumerate([phi_plus2, phi_minus2, psi_plus2, psi_minus2]):
        dm = KetBra(psi_123,psi_123)   
        # Measure
        measure1 = proj(psi2)
        psi_123b = measure1 * psi_123
        
        # Convert collapsed state to density matrix and trace over all qubits except the repeated qubit
        dm = KetBra(psi_123b, psi_123b)
        dm=dm.expand()
        dm = qa.OperatorTrace(dm, over_space=hs[0])
        dm = qa.OperatorTrace(dm, over_space=hs[1])
        dm = qa.OperatorTrace(dm, over_space=hs[2])
        dm = qa.OperatorTrace(dm, over_space=hs[3])
        
        if show_intermediate:
            dm_sympy = convert_to_sympy_matrix(dm.doit())
            triples = dm_sympy.eigenvects()
            for i, triple in enumerate(triples):
                eval1, mult, evecs = triple
                if eval1 != 0:
                    assert len(evecs)==1
                    print(f'{i1} {i2} {evecs[0]}' )
        # Depending on the measurement outcome, the receiver applies a gate to the link qubit.
        if i1==0:
            if i2 == 0:
                gate = None
            elif i2 == 1:
                gate = sigmaz
            elif i2 == 2:
                gate = sigmax
            else:
                gate = sigmay
        elif i1==1:
            if i2 == 0:
                gate = sigmaz
            elif i2 == 1:
                gate = None
            elif i2 == 2:
                gate = sigmay
            else:
                gate = sigmax

        elif i1 == 2:
            if i2 == 0:
                gate = sigmax
            elif i2 == 1:
                gate = sigmay
            elif i2 == 2:
                gate = None
            else:
                gate = sigmaz
        elif i1 == 3:
            if i2 == 0:
                gate = sigmay
            elif i2 == 1:
                gate = sigmax
            elif i2 == 2:
                gate = sigmaz
            else:
                gate = None
        if gate is not None:
            dm = gate(hs[4]).dag() * dm * gate(hs[4])

        # Examine the reconstructed qubit 
        print(f'{i1} {i2} {get_state_from_pure_dm(b*dm,b)}')
        
        # Verify the the reconstructed qubit matches what was transmitted
        x = dm_expected_final/16
        x = x.expand_in_basis()
        y = dm.expand_in_basis()
        assert (x==y)
end=time.time()
print(end-start)            
