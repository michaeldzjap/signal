// lowpass example
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.05, 0.1], [1, 0], [0.01, 0.001]);

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n, f, a, w);
#h, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [h.abs.ampdb, h.phase.unwrap].lace(h.size*2).plot(numChannels: 2);
)

// highpass example
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.02, 0.05], [0, 1], [0.001, 0.01]);

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n, f, a, w);
#h, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [h.abs.ampdb, h.phase.unwrap].lace(h.size*2).plot(numChannels: 2);
)

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

// bandstop example
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.15, 0.3, 0.6, 0.65], [1, 0, 1], [0.01, 0.001, 0.01]);

// fix bump in magnitude response (comment out this line to see the difference)
#f, a, w = FIR.convertToConstraint(f, [1, 0, 1], [0.01, 0.001, 0.01], 0.0005, \type_a);

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n, f, a, w);
#h, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [h.abs.ampdb, h.phase.unwrap].lace(h.size*2).plot(numChannels: 2);
)

// five-band filter
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.17, 0.23, 0.47, 0.53, 0.67, 0.73, 0.82, 0.88], [0, 1, 0, 1, 0], [0.001, 0.01, 0.001, 0.01, 0.001]);

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n + 4, f, a, w, impType: 1);
#h, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [h.abs.ampdb, h.phase.unwrap].lace(h.size*2).plot(numChannels: 2);
)

// another five-band filter
(
var n, f, a, w, h, err, p;

// estimate filter specs
#n, f, a, w = FIR.estimateOrder([0.1, 0.15, 0.3, 0.35, 0.75, 0.8, 0.85, 0.9], [1, 0, 1, 0, 1], [0.01, 0.001, 0.01, 0.001, 0.01]);

// compute FIR filter using Parks/McClellan design algorithm
#h, err = FIR.parksMcClellan(n + 4, f, a, w, impType: 1);
#h, w = h.freqz(n: 300);

// plot magnitude and phase response
p = [h.abs.ampdb, h.phase.unwrap].lace(h.size*2).plot(numChannels: 2);
)