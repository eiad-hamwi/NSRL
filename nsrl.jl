using Beamlines, BeamTracking


Q = 73.0
NN = 181
E0 = 168.5146860368865e9 / NN
EK = 0.3438772962930509e9
EE = (EK + E0) * NN

Ta = Species("181 tantalum +73", Q, E0 * NN, 0.5, 0.0, 0.0, 0.0, HADRON)


R = BeamTracking.E_to_R(Ta, EE)


@eles begin
LD6 = 2.3
BendD6 = -0.15555
BD6 = SBend(L = LD6, angle = BendD6)

L_18D36 = 0.95758
BendD = -0.116355
bend1 = SBend(L = L_18D36, angle = BendD, e1 = BendD/2, e2 = BendD/2)
bend2 = SBend(L = L_18D36, angle = BendD, e1 = BendD/2, e2 = BendD/2)
bend3 = SBend(L = L_18D36, angle = BendD, e1 = BendD/2, e2 = BendD/2)

LNSRLQh = 0.3555

Q1 = Quadrupole(L = LNSRLQh * 2, Bn1 = -2.32)
Q2 = Quadrupole(L = LNSRLQh * 2, Bn1 = 2.7)
Q3 = Quadrupole(L = LNSRLQh * 2, Bn1 = 2.95)
Q4 = Quadrupole(L = LNSRLQh * 2, Bn1 = -2.97)
Q5 = Quadrupole(L = LNSRLQh * 2, Bn1 = 1.87)
Q6 = Drift(L = LNSRLQh * 2)
Q7 = Quadrupole(L = LNSRLQh * 2, Bn1 = -0.99)
Q8 = Quadrupole(L = LNSRLQh * 2, Bn1 = 1.08)    

L_oct = 0.6096

O1 = Octupole(L = L_oct, Bn3 = -1.75E+03)
O2 = Octupole(L = L_oct, Bn3 = 1.75E+03)


D1 = Drift(L = 15.454)
D2 = Drift(L = 1.873)
D3 = Drift(L = 10.335)
D4 = Drift(L = 0.363)
D5 = Drift(L = 0.312)
D6 = Drift(L = 3.212)
D7 = Drift(L = 1.873)
D8 = Drift(L = 11.003)
D9 = Drift(L = 0.328)
D10 = Drift(L = 7.091)
D11 = Drift(L = 0.326)
D12 = Drift(L = 6.173)
D13 = Drift(L = 5.94)
D14 = Drift(L = 24.19)
end

nsrl = Beamline([BD6, D1, Q1, D2, Q2, D3, bend1, D4, bend2, D5, bend3, D6, Q3, 
D7, Q4, D8, O1, D9, Q5, D10, O2, D11, Q6, D12, Q7, D13, Q8, D14], 
R_ref = R, species_ref = Ta)


for ele in nsrl.line
    ele.aperture_shape = ApertureShape.Rectangular
    ele.x1_limit = -0.2
    ele.x2_limit = 0.2
    ele.y1_limit = -0.2
    ele.y2_limit = 0.2
end