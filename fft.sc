+ ArrayedCollection {

	/**
	 * Chirp z-transform
	 * This function returns the m-element z-transform of sequence x,
     * where m, w and a are scalars which specify the contour in the z-plane
     * on which the z-transform is computed. m is the length of the transform,
     * w is the complex ratio between points on the contour, and a is the
     * complex starting point. More explicitly, the contour in the z-plane
     * (a spiral or "chirp" contour) is described by
	 *     z = a * (w**(-(0..m - 1)))
     *
     * The parameters m, w, and a are optional; their default values are
	 * m = length(x), w = exp(Complex(0, 1)*2pi/m), and a = 1. These defaults
     * cause czt to return the z-transform of x at equally spaced points
     * around the unit circle, equivalent to x.fft.
	 *
	 * NOTE: this is a SC port of MATLAB's czt function
	 */
	czt { arg k = this.size, w, a = 1;
		var nfft, m, kk, kk2, ww, nn, aa, x, y, fv;

		m = this.size;
		w = w ? exp(Complex(0, -1)*2pi/k);
		nfft = (m + k - 1).nextPowerOfTwo;   // length for power of two fft

		// pre-multiply data
		kk = (m.neg + 1,m.neg + 2..(k - 1).max(m - 1));
		kk2 = kk**2*0.5;
		ww = kk2 collect: { |p| w**p };
		nn = (0..m - 1);
		aa = (nn.size collect: { |i|
			var wi = ww[m - 1 + nn[i]];
			a**nn[i].neg*[wi.real, wi.imag]
		});
		y = (this*aa).flat;

		// fast convolution via FFT
		y = y.zeropad(nfft*2);
		y = y.cfft(nfft, true);
		fv = ((k - 1 + m) collect: { |i| var wi = 1/ww[i]; [wi.real, wi.imag] }).flat;
		fv = fv.zeropad(nfft*2);
		fv = fv.cfft(nfft, true);   // Chirp filter
		y = (y.size.div(2) collect: { |i| i = i*2;   // complex multiplication
			[y[i]*fv[i] - (y[i + 1]*fv[i + 1]), y[i]*fv[i + 1] + (y[i + 1]*fv[i])]
		}).flat;
		y = y.cfft(nfft, true, false)*nfft*2;
		y = y.clump(2) collect: { |c| Complex(c[0], c[1]) };

		// final multiply
		^(y[m - 1..m + k - 2]*ww[m - 1..m + k - 2]);
	}

	/**
	 * rfft replaces 2*n real data points in x with n complex
	 * values representing the positive frequency half of their
	 * Fourier spectrum. n must be a power of 2.
	 *
	 * source: Moore, F. Richard. Elements of Computer Music,
	 * P T R Prentice Hall, 1990.
	 */
	rfft { arg n, inplaceFlag = false, forward = true;
		var c1, c2, theta, wr, wi, x, xr, xi, wpr, wpi, temp, n2p1, i, i1, i2, i3, i4, h1r, h1i, h2r, h2i;

		x = inplaceFlag.if { this } { this.copy };

		theta = pi*n.reciprocal;
		wr = 1;
		wi = 0;
		c1 = 0.5;
		forward.if {
			c2 = -0.5;
			x.cfft(n, true);
			xr = x[0];
			xi = x[1];
		} {
			c2 = 0.5;
			theta = theta.neg;
			xr = x[1];
			xi = 0;
			x[1] = 0;
		};
		wpr = -2*(sin(0.5*theta)**2);
		wpi = sin(theta);
		n2p1 = (n << 1) + 1;
		i = 0;
		while({ i <= (n >> 1) }) {
			i1 = i << 1;
			i2 = i1 + 1;
			i3 = n2p1 - i2;
			i4 = i3 + 1;
			(i == 0).if {
				h1r = c1*(x[i1] + xr);
				h1i = c1*(x[i2] - xi);
				h2r = c2.neg*(x[i2] + xi);
				h2i = c2*(x[i1] - xr);
				x[i1] = h1r + (wr*h2r) - (wi*h2i);
				x[i2] = h1i + (wr*h2i) + (wi*h2r);
				xr = h1r - (wr*h2r) + (wi*h2i);
				xi = h1i.neg + (wr*h2i) + (wi*h2r);
			} {
				h1r = c1*(x[i1] + x[i3]);
				h1i = c1*(x[i2] - x[i4]);
				h2r = c2.neg*(x[i2] + x[i4]);
				h2i = c2*(x[i1] - x[i3]);
				x[i1] = h1r + (wr*h2r) - (wi*h2i);
				x[i2] = h1i + (wr*h2i) + (wi*h2r);
				x[i3] = h1r - (wr*h2r) + (wi*h2i);
				x[i4] = h1i.neg + (wr*h2i) + (wi*h2r);
			};
			temp = wr;
			wr = temp*wpr - (wi*wpi) + wr;
			wi = wi*wpr + (temp*wpi) + wi;

			i = i + 1;
		};
		x[1] = xr; // this replaces x[1] with the real part of the Nyquist frequency value

		^x.sanitise
	}

	/**
	 * cfft replaces/returns array x containing nc complex values
	 * (2*nc values alternating real, imaginary, and so on)
	 * by its Fourier transform.
	 *
	 * source: Moore, F. Richard. Elements of Computer Music,
	 * P T R Prentice Hall, 1990.
	 */
	cfft { arg nc, inplaceFlag = false, forward = true;
		var wr, wi, wpr, wpi, theta, mmax, nd, m, i, j, delta, x;

		x = inplaceFlag.if { this } { this.copy };
		nd = nc << 1;

		x.bitreverse(nd);

		mmax = 2;
		while({ mmax < nd }) {
			delta = mmax << 1;
			theta = 2pi/(forward.if { mmax } { mmax.neg });
			wpr = -2*(sin(0.5*theta)**2);
			wpi = sin(theta);
			wr = 1;
			wi = 0;
			m = 0;
			while({ m < mmax }) {
				var rtemp, itemp;

				i = m;
				while({ i < nd }) {
					j = i + mmax;
					rtemp = wr*x[j] - (wi*x[j + 1]);
					itemp = wr*x[j + 1] + (wi*x[j]);
					x[j] = x[i] - rtemp;
					x[j + 1] = x[i + 1] - itemp;
					x[i] = x[i] + rtemp;
					x[i + 1] = x[i + 1] + itemp;

					i = i + delta;
				};
				rtemp = wr;
				wr = rtemp*wpr - (wi*wpi) + wr;
				wi = wi*wpr + (rtemp*wpi) + wi;

				m = m + 2;
			};

			mmax = delta;
		};

		// scale output
		x = x*(forward.if { nd.reciprocal } { 2 });

		^x.sanitise
	}

	/**
	 * bitreverse places 'this' containing n/2 complex values
	 * into bit-reversed order
	 *
	 * source: Moore, F. Richard. Elements of Computer Music,
	 * P T R Prentice Hall, 1990.
	 */
	bitreverse { arg n;
		var rtemp, itemp, i, j, m;

		i = 0; j = 0;
		while({ i < n }) {
			(j > i).if {
				rtemp = this[j]; itemp = this[j + 1];   // complex exchange
				this[j] = this[i]; this[j + 1] = this[i + 1];
				this[i] = rtemp; this[i + 1] = itemp;
			};

			m = n >> 1;
			while({ m >= 2 and: { j >= m } }) { j = j - m; m = m >> 1; };

			i = i + 2; j = j + m;
		}
	}

	// same as bitreverse, but for real values
	bitreverse2 {
		var temp, i, j, m, n;

		n = this.size;

		i = 0; j = 0;
		while({ i < n }) {
			(j > i).if {
				temp = this[j];
				this[j] = this[i];
				this[i] = temp;
			};

			m = n >> 1;
			while({ m >= 2 and: { j >= m } }) { j = j - m; m = m >> 1; };

			i = i + 1; j = j + m;
		}
	}

	sanitise { arg eps = 1e-9;
		this.size do: { |i|
			(this[i].abs < 1e-9).if { this[i] = 0.0 }
		};
		^this
	}

	// zero pad/truncate
	zeropad { arg n;
		var l = this.size, x = this;
		(l < n).if {   // zero-pad data
			x = x ++ (0!(n - l));
		};
		(l > n).if {   // truncate data
			(l - n) do: { |i| x.removeAt(l - 1 - i) }
		};
		^x
	}

}