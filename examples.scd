// bandpass example
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.2, 0.25, 0.6, 0.7], [0, 1, 0], [0.001, 0.01, 0.01]);

// fix bump in magnitude response (comment out this line to see the difference)
#f, a, w = FIR.convertToConstraint(f, [0, 1, 0], [0.001, 0.01, 0.01], 0.0005, \type_a);

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n, f, a, w);
#h, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [h.abs.ampdb, h.phase.unwrap].lace(h.size*2).plot(numChannels: 2);
)