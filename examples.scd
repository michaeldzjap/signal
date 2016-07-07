// bandpass example
(
var n, f, a, w, h, err;
#n, f, a, w = FIR.estimateOrder([0.2, 0.25, 0.6, 0.7], [0, 1, 0], [0.001, 0.01, 0.01]);
//#f, a, w = FIR.convertToConstraint(f, [0, 1, 0], [0.001, 0.01, 0.01], 0.0005, \type_a);
#h, err = FIR.parksMcClellan(n, f, a, w);
h.freqz(n: 20);
)

h = [[1, 0.5], [1, -0.25]].freqz(n: 20)

[h[0].abs, h[0].phase].lace(h[0].size*2).plot(numChannels: 2)

[[1, 2], [3, 4]].lace(4)