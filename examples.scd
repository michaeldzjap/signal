// bandpass example
(
var n, f, a, w, h, err;
#n, f, a, w = FIR.estimateOrder([0.2, 0.25, 0.6, 0.7], [0, 1, 0], [0.001, 0.01, 0.01]);
#f, a, w = FIR.convertToConstraint(f, [0, 1, 0], [0.001, 0.01, 0.01], 0.0005, \type_a);
#h, err = FIR.parksMcClellan(n, f, a, w);
#h, w = h.freqz(n: 512);
[h.abs, h.phase].lace(h.size*2).plot(numChannels: 2)
)