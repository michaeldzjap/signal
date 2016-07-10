/**
 * Code for designing optimum FIR filters using the Parks-McClellan algorithm.
 *
 * This code is more or less a straight translation of the MATLAB function firpm to
 * SuperCollider, with the exception that the Remez exchange part of the algorithm is
 * replaced by the implementation(s) proposed in [1].
 *
 * [1] Ahsan, M. & Saramäki, T. (2012). Two Novel Implementations of the Remez Multiple
 *     Exchange Algorithm for Optimum FIR Filter Design, In:MATLAB – A Fundamental Tool
 *     for Scientific Computing and Engineering Applications, Vol. 2, Katsikis, V.N., (Eds.),
 *     Ch. 3, (Sep. 2012), InTech, pp. 37-66.
 */

FIR {

	*parksMcClellan { arg order, freqs, amps, wtx, ftype, lgrid, impType = 0;
		var nfilt, grid, des, wt, signVal, hilbert, neg, h, err, iext;

		#nfilt, freqs, grid, des, wt, ftype, signVal, hilbert, neg = FIR.prInit(order, freqs, amps, wtx, ftype, lgrid);

		#h, err, iext = FIR.prRemez(nfilt, freqs/2, grid/2, des, wt, neg, impType);

		err = err.abs;
		h = h ++ ((0.5 - neg).sign*h[h.lastIndex - (nfilt%2)..0]);
		h = h[h.lastIndex..0];
		(neg.booleanValue and: { hilbert.not }).if { h = h.neg }; // make sure differentiator has correct sign

		^[h, err]
	}

	*prInit { arg order, freqs, amps, wtx, ftype, lgrid;
		var nargs = 0, df, exceptionFlag, msg1, msg2, signVal, nfilt, nodd, hilbert, neg, grid, des, wt;

		// Check nr. of input arguments
		[order, freqs, amps, wtx, ftype, lgrid].checkNumArgs(3, 6);

		(order < 3).if { Error("Filter order must be 3 or more.").throw };

		/**
		 * Define default values for optional input arguments:
		 */
		wtx = wtx ? (1 ! floor((1 + freqs.size)/2));
		ftype = ftype ? \f;
		lgrid = lgrid ? 16;

		/**
		 * Error checking
		 */
		(freqs.size % 2 != 0).if { Error("The number of frequency points must be even.").throw };
		(freqs any: { |f| f < 0 } or: { freqs any: { |f| f > 1 } }).if {
			Error("Normalised frequencies must lie between 0 and 1.").throw
		};
		df = freqs.differentiate; df.removeAt(0);
		(df any: { |f| f < 0 }).if { Error("Frequencies must be non-decreasing.").throw };
		(wtx.size != (1 + freqs.size).div(2)).if { Error("There should be one weight per band.").throw };
		((wtx any: { |w| w.sign == 1 } and: { wtx any: { |w| w.sign == -1 } })
			or: { wtx any: { |w| w.sign == 0 } }).if {
			Error("All weights must be greater than zero.").throw
		};

		/**
		 * Check for valid filter length
		 */
		exceptionFlag = 0;
		#[\differentiator, \hilbert, \h, \d].includes(ftype).if { exceptionFlag = 1 };
		#order, msg1, msg2 = FIR.prCheckOrder(order, freqs.last, amps.last, exceptionFlag);
		msg1.notNil.if { Error(msg1).throw };
		msg2.notNil.if { (msg2 ++ "Alternatively, you can pass a trailing \\h argument, as in firpm(N, F, A, W, \\h), to design a type 4 linear phase filter.").warn };

		/**
		 * Determine symmetry of filter
		 */
		signVal = 1;
		nfilt = order + 1;      // filter length
		nodd = nfilt % 2;       // nodd == 1 ==> filter length is odd, nodd == 0 ==> filter length is even

		hilbert = false;
		(ftype.toLower == \h or: { ftype.toLower == \hilbert }).if {
			ftype = 3;    // Hilbert transformer
			hilbert = true;
			nodd.booleanValue.not.if { ftype = 4 };
		};
		(ftype.toLower == \d or: { ftype.toLower == \differentiator }).if {
			ftype = 4;    // Differentiator
			signVal = -1;
			nodd.booleanValue.if { ftype = 3 };
		} {
			ftype = 1;                  // Regular filter
			nodd.booleanValue.not.if { ftype = 2 };
		};

		// neg == 1 ==> antisymmetric imp resp, neg == 0 ==> symmetric imp resp
		(ftype == 3 or: { ftype == 4 }).if { neg = 1 } { neg = 0 };

		/**
		 * Create grid of frequencies on which to perform Remez exchange iteration
		 */
		grid = FIR.createFrequencyGrid(nfilt, lgrid, freqs, neg, nodd);
		while({ grid.size <= nfilt }) {
			lgrid = lgrid*4;   // need more grid points
			grid = FIR.createFrequencyGrid(nfilt, lgrid, freqs, neg, nodd);
		};

		/**
		 * Get desired frequency characteristics at the frequency points
		 * in the specified frequency band intervals.
		 */
		#des, wt = FIR.prFrequencyResponse(order, freqs, grid, wtx, amps, ftype == 4);

		^[nfilt, freqs, grid, des, wt, ftype, signVal, hilbert, neg]
	}

	// Check if specified filter order is valid.
	*prCheckOrder { arg order, freqEnd, ampEnd, exceptionFlag;
		var msg1, msg2, oddord = false;

		[order, freqEnd, ampEnd, exceptionFlag].checkNumArgs(3, 4);

		exceptionFlag.isNil.if { exceptionFlag = false };

		(order.isNil or: { order.isInteger.not } or: { order <= 0 }).if {
			msg1 = "Filter order must be a real, positive integer.";
			^[order, msg1, msg2]
		};

		(order % 2 == 1).if { oddord = true };

		(ampEnd != 0 and: { freqEnd == 1 } and: { oddord } and: { exceptionFlag.not }).if {
			msg2 = "Odd order symmetric FIR filters must have a gain of zero at the Nyquist frequency. The order is being increased by one.";
			order = order + 1;
		};

		^[order, msg1, msg2]
	}

	*prFrequencyResponse { arg order, freqs, gridFreqs, weights, amps, diffFlag = false;
		var dh = [], dw = [], nbands, l;

		// Prevent discontinuities in desired function
		(1, 3 .. freqs.size - 3) do: { |k|
			(freqs[k] == freqs[k + 1]).if { Error("Adjacent bands not allowed.").throw };
		};
		(freqs.size != amps.size).if { Error("Frequency and amplitude arrays must be the same length.").throw };

		nbands = amps.size/2;
		l = 0;
		while({ (l + 1)/2 <= nbands }) {
			var selFreqs = gridFreqs select: { |gf| gf >= freqs[l] and: { gf <= freqs[l + 1] } };
			// desired magnitude is line connecting amps[l] to amps[l + 1]
			(freqs[l + 1] != freqs[l]).if {
				var slope = (amps[l + 1] - amps[l])/(freqs[l + 1] - freqs[l]);
				dh = dh ++ (selFreqs collect: { |sf| slope*sf + amps[l] - (slope*freqs[l]) });
			} {   // zero bandwidth band
				dh = dh ++ ({ (amps[l] + amps[l + 1])/2 } ! selFreqs.size);
			};
			dw = dw ++ (weights[(l + 1).div(2)]/(1 + ((diffFlag and: { amps[l + 1] >= 1e-4 }).asInteger*(selFreqs/2 - 1))));
			l = l + 2;
		};
		^[dh, dw]
	}

	*createFrequencyGrid { arg nfilt, lgrid, freqs, neg, nodd;
		var grid = List(), nfcns, delf, eps = 2.220446049250313e-016, j, l, ngrid;
		nfcns = nfilt.div(2);
		(nodd == 1 and: { neg == 0 }).if { nfcns = nfcns + 1 };
		grid.add(freqs.first);
		delf = (lgrid*nfcns).reciprocal;
		// If value at frequency 0 is contrained, make sure first grid point is not too small
		(neg != 0 and: { grid.first < delf }).if {
			// Handle narrow bands
			(freqs.first > sqrt(eps)).if {
				grid[0] = freqs.first
			} {
				(delf < freqs[1]).if { grid[0] = delf } { grid[0] = 0.5*(freqs[1] - freqs.first) }
			}
		};
		j = 0; l = 0;
		while({ l + 1 <= freqs.lastIndex }) {
			var fup, newgrid, jend;
			fup = freqs[l + 1];
			newgrid = (grid[j] + delf, grid[j] + delf + delf .. fup + delf);
			(newgrid.size < 11).if {
				var delf1 = (fup + delf - (grid[j] + delf))/10;
				newgrid = (grid[j] + delf1, grid[j] + delf1 + delf1 .. fup + delf1);
			};
			grid = grid ++ newgrid;
			jend = grid.lastIndex;
			(jend > 0).if { grid[jend - 1] = fup; j = jend } { j = jend + 1 };

			l = l + 2;
			(l + 1 <= freqs.lastIndex).if { grid[j] = freqs[l] };
		};
		ngrid = j - 1;
		// If value at frequency 1 is contrained, remove that grid point
		(neg == nodd and: { grid[ngrid] > (1 - delf) }).if {
			(freqs[freqs.lastIndex - 1] < (1 - delf)).if {
				ngrid = ngrid - 1;
			} {
				grid[ngrid] = freqs[freqs.lastIndex - 1];
			}
		};
		grid = grid[0..ngrid];
		^grid
	}

	*prRemez { arg nfilt, edge, grid, des, wt, neg, impType;
		var nbands, jb, nodd, nfcns, ngrid, temp, j, iext, nz, x, y, ad, dev, nm1, fsh, cn, delf, l, kkk, a, aa, bb, dden, alpha, h;

		nbands = edge.size/2;
		jb = 2*nbands;
		nodd = nfilt % 2;
		nfcns = nfilt.div(2);
		(nodd == 1 and: { neg == 0 }).if { nfcns = nfcns + 1 };

		ngrid = grid.size;

		(neg <= 0).if {
			(nodd != 1).if {
				des = des/cos(pi*grid);
				wt = wt*cos(pi*grid);
			};
		} {
			(nodd != 1).if {
				des = des/sin(pi*grid);
				wt = wt*sin(pi*grid);
			} {
				des = des/sin(2*pi*grid);
				wt = wt*sin(2*pi*grid);
			}
		};
		temp = (ngrid - 1)/nfcns;
		j = (0..nfcns - 1);
		iext = floor(temp*j ++ [ngrid - 1]);
		nz = nfcns + 1;

		#x, y, ad, dev = (impType == 0).if { FIR.prRemezImp1(nz, iext, ngrid, grid, des, wt) } { FIR.prRemezImp2(nz, iext, ngrid, grid, des, wt) };

		// Inverse Fourier transformation
		nm1 = nfcns - 1;
		fsh = 1e-6;
		x = x ++ [-2];
		cn = 2*nfcns - 1;
		delf = cn.reciprocal;
		l = 0;
		kkk = false;
		(edge.first == 0 and: { edge[jb] == 0.5 } or: { nfcns <= 3 }).if {
			kkk = true;
		};
		kkk.not.if {
			var dtemp, dnum;
			dtemp = cos(2pi*grid.first);
			dnum = cos(2pi*grid[ngrid - 1]);
			aa = 2/(dtemp - dnum);
			bb = (dtemp + dnum).neg/(dtemp - dnum);
		};
		a = 0!nfcns;
		nfcns do: { |j|
			var ft, xt, xe;
			ft = j*delf;
			xt = cos(2pi*ft);
			kkk.not.if {
				xt = (xt - bb)/aa;
				ft = acos(xt)/2pi;
			};
			xe = x[l];
			while({ xt <= xe and: { xe - xt >= fsh } }) {
				l = l + 1;
				xe = x[l];
			};
			((xt - xe).abs < fsh).if {
				a[j] = y[l]
			} {
				var c;
				grid[0] = ft;
				// gee
				c = ad/(cos(2pi*ft) - x[0..nz - 1]);
				a[j] = (c*y).sum/c.sum;
			};
			l = 0.max(l - 1);
		};

		dden = 2pi/cn;

		alpha = 0!nfcns;
		nfcns do: { |j|
			var dnum;
			dnum = j*dden;
			(nm1 < 1).if {
				alpha[j] = a[0];
			} {
				alpha[j] = a[0] + (2*cos(dnum*(1..nm1))*a[1..nfcns - 1]).sum;
			}
		};
		alpha = ([alpha.first] ++ (2*alpha[1..nfcns - 1])/cn);
		kkk.not.if {
			var p = 0!2, q = 0!1;
			p[0] = 2*alpha[nfcns - 1]*bb + alpha[nm1 - 1];
			p[1] = 2*aa*alpha[nfcns - 1];
			q[0] = alpha[nfcns - 3] - alpha[nfcns - 1];
			(nm1 - 1) do: { |j| var sel, jm1, jp1; j = j + 1;
				(j == (nm1 - 1)).if {
					aa = aa/2;
					bb = bb/2;
				};
				p = p ++ [0];
				sel = (0..j);
				a.overWrite(p|@|sel, sel.first);
				p.overWrite(2*bb*(a|@|sel), sel.first);
				p[1] = p[1] + (2*a.first*aa);
				jm1 = j - 1;
				sel = (0..jm1);
				p.overWrite((p|@|sel) + (q|@|sel) + (aa*(a|@|(sel + 1))), sel.first);
				jp1 = j + 1;
				sel = (2..jp1);
				p.overWrite((p|@|sel) + (aa*(a|@|(sel - 1))), sel.first);
				(j != (nm1 - 1)).if {
					sel = (0..j);
					q = q ++ [0];
					q.overWrite((a|@|sel).neg, sel.first);
					q[0] = q.first + alpha[nfcns - 3 - j];
				}
			};
			alpha.overWrite(p[0..nfcns - 1], 0);
		};
		(nfcns <= 3).if {
			alpha[nfcns] = 0;
			alpha[nfcns + 1] = 0;
		};

		(neg <= 0).if {
			(nodd != 0).if {
				h = 0.5*alpha[nz - 2..nz - nm1 - 1] ++ [alpha[0]];
			} {
				h = 0.25*([alpha[nfcns - 1]] ++ (alpha[nz - 3..nz - nm1 - 1] + alpha[nfcns - 1..nfcns - nm1 + 1]) ++ [2*alpha[0] + alpha[1]]);
			}
		} {
			(nodd != 0).if {
				h = 0.25*([alpha[nfcns - 1]] ++ [alpha[nm1 - 1]]);
				h = h ++ (0.25*((alpha[nz - 4..nz - nm1 - 1] - alpha[nfcns - 1..nfcns - nm1 + 2]) ++ [2*alpha[0] - alpha[2]])) ++ [0];
			} {
				h = 0.25*([alpha[nfcns - 1]] ++ (alpha[nz - 3..nz - nm1 - 1] - alpha[nfcns - 1..nfcns - nm1 + 1]) ++ [2*alpha[0] - alpha[1]]);
			}
		};

		^[h, dev, iext]
	}

	*prGenWeightedError { arg nz, ngrid, grid, des, wt, l_trial;
		var x, y, a, ad = 1!nz, dev, add = 1!nz, dnum = 0, dden = 0, x_all, err_num, err_den, err_cy, wei_err, dev_vect;

		/**
		 * Step I: based on the present 'trial' vector l_trial, generate the
		 * weighted error function wei_err(k) at all the grid points.
		 */
		x = cos(2pi*(grid|@|l_trial));   // step 1: Lagrange abscissa vector x
		a = Matrix.withFlatArray(x.size, 1, x)*Matrix.withFlatArray(1, nz, 1!nz) - (Matrix.withFlatArray(nz, 1, 1!nz)*Matrix.withFlatArray(1, x.size, x));
		nz do: { |k|
			a[k, k] = 1.0;
			a.doCol(k, { |entree| ad[k] = ad[k]*entree });
		};
		ad = ad*((-2)**(nz - 1));  // step 1: Lagrange coefficient vector ad...
		ad = ad.reciprocal;        // found efficiently without using the function remezdd
		nz.div(2) do: { |i| add[i*2 + 1] = -1 };
		ad.size do: { |i| dnum = dnum + (ad[i]*(des|@|l_trial)[i]) };
		add.size do: { |i| dden = dden + (add[i]*(ad/(wt|@|l_trial))[i]) };
		dev = dnum.neg/dden;   // step 1: current value of deviation
		// step 2: Lagrange ordinate vector y
		y = des|@|l_trial + (dev*add/(wt|@|l_trial));
		// step 3: overall obscissa vector x_all
		x_all = cos(2pi*grid[0..ngrid - 1]);
		err_num = 0!ngrid;    // step 4: initialisaztion of err_num + err_den
		err_den = 0!ngrid;
		nz do: { |jj|         // step 5 and 6: intermediate evaluations for...
			var aid = ad[jj]/(x_all - x[jj]);   // obtaining the weighted error...
			err_den = err_den + aid;            // wei_err[k] at all the grid points.
			err_num = err_num + (y[jj]*aid);
		};
		err_cy = err_num/err_den;
		wei_err = (err_cy - des)*wt;                      // step 7: generate the vector wei_err
		dev_vect = 1!l_trial.size;
		l_trial.size.div(2) do: { |i| dev_vect[i*2 + 1] = dev_vect[i*2 + 1].neg };
		dev_vect = dev_vect*dev;                           // entries of wei_err at l_trial[0:nz - 1]...
		l_trial do: { |lt, i| wei_err[lt] = dev_vect[i] }; // by using the values of dev (-dev).
		^[wei_err, dev_vect, l_trial, x, y, ad, dev]
	}

	*prRemezImp1 { arg nz, iext, ngrid, grid, des, wt;
		var niter, itrmax, l_trial, a, x, y, ad, dev;

		// Initialisation phase
		niter = 1;                  // Initialise the iteration counter
		itrmax = 250;               // Maximum number of iterations
		l_trial = iext[0..nz - 1];  // Startup value of l_trial

		// Iteration phase
		// Remez loop for locating desired nz indices among the grid points
		block { |break|
			while({ niter < itrmax }) {
				var /*add = 1!nz, dnum = 0, dden = 0, x_all, err_num, err_den, err_cy,*/ wei_err, dev_vect, l_real=nil!nz, err_vic = 0!nz, endsearch, err_end = 0!nz, l_end_real = 0!nz;
				ad = 1!nz;

				#wei_err, dev_vect, l_trial, x, y, ad, dev = FIR.prGenWeightedError(nz, ngrid, grid, des, wt, l_trial);

				// Step II: perform vicinity search
				nz do: { |k|   // steps 1, 2 and 3: start of vicinity search
					var ind_vicinity, err_vicinity, low;
					(k == 0).if {
						low = 0;
						err_vicinity = (wei_err[l_trial[k]]).sign*(wei_err|@|(low..l_trial[1] - 1));
					} {
						(k == (nz - 1)).if {
							low = (l_trial[k - 1] + 1).max(l_real[k - 1] + 1);
							err_vicinity = wei_err[l_trial[k]].sign*(wei_err|@|(low..ngrid - 1));
						} {
							low = (l_trial[k - 1] + 1).max(l_real[k - 1] + 1);
							err_vicinity = wei_err[l_trial[k]].sign*(wei_err|@|(low..l_trial[k + 1] - 1));
						}
					};
					ind_vicinity = err_vicinity.maxIndex;

					l_real[k] = ind_vicinity + low;
					(k == 0 or: { k == (nz - 1) }).if {   // step 3: find err_vic[0] = wei_err[l_real[0]]...
						err_vic[k] = wei_err[l_real[k]]
					}
				};

				// Step III: perform endpoint search
				endsearch = 0;    // step 1: start endpoint search
				err_end[0] = 0;   // step 2: needed for the case, where upp_end = 0
				(l_real[0] > 1 and: { l_trial[0] > 1 }).if {   // step 2: find l_end_true[0]...
					var upp_end, err_endpoint, ind_endpoint;
					upp_end = (l_real[0] - 1).min(l_trial[0] - 1);   // and err_end[0]
					err_endpoint = (wei_err[l_real[0]]).sign.neg*(wei_err|@|(0..upp_end));
					ind_endpoint = err_endpoint.maxIndex;
					l_end_real[0] = ind_endpoint;
					err_end[0] = (wei_err[l_real[0]]).sign.neg*wei_err[l_end_real[0]];
					(err_end[0] > err_vic[nz - 1].abs).if { endsearch = 1; }; // step 3: use 'endsearch=1' or not?
				};
				(l_real[nz - 1] < (ngrid - 1) and: { l_trial[nz - 1] < (ngrid - 1) }).if {   // step 4: find...
					var low_end, err_endpoint, ind_endpoint;
					low_end = (l_real[nz - 1] + 1).max(l_trial[nz - 1] + 1);   // l_end_real[nz - 1]...
					err_endpoint = (wei_err[l_real[nz - 1]]).sign.neg*(wei_err|@|(low_end..ngrid - 1));
					ind_endpoint = err_endpoint.maxIndex;
					l_end_real[nz - 1] = ind_endpoint + low_end;
					err_end[nz - 1] = (wei_err[l_real[nz - 1]]).sign.neg*wei_err[l_end_real[nz - 1]];
					(err_end[nz - 1] > err_vic[0].abs.max(err_end[0])).if {   // step 5:...
						endsearch = 2;   // use 'endsearch=2' or not?
					}
				};
				(endsearch == 1).if {   // step 7: 'endsearch=1' is valid. Form...
					l_real = [l_end_real[0]] ++ l_real[0..nz - 2];   // l_real accordingly
				};
				(endsearch == 2).if {   // step 8: 'endsearch=2' is true. Form...
					l_real = l_real[1..nz - 1] ++ [l_end_real[nz - 1]];   // l_real accordingly
				};

				// Segment 4: test convergence
				(l_real == l_trial).if {   // step 1: the real and trial vectors coincide.
					break.value;           // Hence, stop. Remez loop ended successfully.
				} {
					l_trial = l_real;      // step 2: otherwise, replace the values of l_trial
					niter = niter + 1;     // with the values of l_real and continue
				}
			}
		};
		^[x, y, ad, dev]
	}

	*prRemezImp2 { arg nz, iext, ngrid, grid, des, wt;
		var niter, itrmax, l_trial, a, x, y, ad, dev;

		// Initialisation phase
		niter = 1;                  // Initialise the iteration counter
		itrmax = 250;               // Maximum number of iterations
		l_trial = iext[0..nz - 1];  // Startup value of l_trial

		// Iteration phase
		// Remez loop for locating desired nz indices among the grid points
		block{ |break|
			while({ niter < itrmax }) {
				var wei_err, dev_vect, l_aid1, l_aid2, ind, rowInds, colInds, weiErrAbs, l_real_start, l_real_init, wei_real, wei_comp, l_real;

				#wei_err, dev_vect, l_trial, x, y, ad, dev = FIR.prGenWeightedError(nz, ngrid, grid, des, wt, l_trial);

				// Step II: determine the vector l_real_start
				// step 1: find l_aid1
				l_aid1 = (wei_err ++ [0]).differentiate.sign.differentiate;
				l_aid1.removeAt(0);
				l_aid1 = l_aid1 selectIndices: { |item| item != 0 };
				// step 2: determine l_aid2
				l_aid2 = l_aid1|@|((wei_err|@|l_aid1).abs selectIndices: { |item| item >= dev.abs });
				// step 3: this is not ideal, but no sparse matrix support in SC
				rowInds = (0..l_aid2.lastIndex);
				colInds = ([1] ++ ([(wei_err|@|(l_aid2[1..l_aid2.lastIndex])) >= 0,(wei_err|@|(l_aid2[0..l_aid2.lastIndex - 1])) >= 0].flop collect: { |item| item[0] != item[1] }).asInteger).integrate - 1;
				ind = Matrix.fill(rowInds.maxItem + 1, colInds.maxItem + 1, { 0 });
				weiErrAbs = (wei_err|@|l_aid2).abs;
				l_aid2.size do: { |i| ind.put(rowInds[i], colInds[i], weiErrAbs[i]) };
				ind = ind.cols collect: { |i| ind.getCol(i).maxIndex }; // retrieve index of max item per col
				// step 4: determine l_real_start
				l_real_start = l_aid2|@|ind;

				// Step III: determine the vector l_real
				l_real_init = l_real_start;   // step 1
				((l_real_init.size - nz)%2 == 1).if {   // step 2: odd difference
					(wei_err[l_real_init.first].abs <= wei_err[l_real_init.last].abs).if {
						l_real_init.removeAt(0);   // step 3: discard the first entry...
					} {   // of l_real_init
						l_real_init.removeAt(l_real_init.lastIndex);   // otherwise discard the last entry
					}
				};
				while({ l_real_init.size > nz }) {   // step 4
					wei_real = (wei_err|@|l_real_init).abs;    // start of step 5
					wei_comp = (wei_real.size - 1) collect: { |i| wei_real[i].max(wei_real[i + 1]) }; // end of step 5
					(wei_err[l_real_init.first].abs.max(wei_err[l_real_init.last].abs) <= wei_comp.minItem).if {  // start of step 6
						l_real_init = l_real_init[1..l_real_init.lastIndex - 1];   // end of step 6
					} {
						var ind_omit = wei_comp.minIndex;   // start: step 7
						2 do: { l_real_init.removeAt(ind_omit) };
					}
				};
				l_real = l_real_init;

				// Step IV: test convergence
				(l_real every: { |lr, i| lr == l_trial[i] }).if {
					break.value;  // step 1: The real and trial vectors coincide. Hence stop. Remez loop ended successfully
				} {
					l_trial = l_real;   // step 2: Otherwise, replace the values of l_trial with the values of l_real and continue
					niter = niter + 1;
				}
			}
		}
		^[x, y, ad, dev]
	}

	/**
	 * This function converts the original design specifications into Type A
	 * or Type B transition band constraints compatible specifications.
	 *
	 * Input parameters:
	 * - fo contains the edges of the r bands
	 * - deso contains the desired values in the r bands
	 * - devo contains the admissable deviations in the r bands
	 * - alpha is a small positive constant
	 * - type=1 (2) generates Type A (B) transition band constraints
	 *
	 * Output parameters:
	 * - f contains the edges of the 2r - 1 bands
	 * - des contains the desired values on all the edges of 2r - 1 bands
	 * - wt contains the weights in the 2r - 1 bands
	 */
	*convertToConstraint { arg fo, deso, devo, alpha, itype;
		var r, f, des, dev, wt;

		// Check if the input data is correct
		(alpha < 0).if { Error("alpha should be a small positive number.").throw };
		r = deso.size;

		// Generate the output parameters
		f = 0!(4*(r - 1) + 2);
		des = 0!(2*(r - 1) + 1);
		dev = 0!(2*(r - 1) + 1);
		r do: { |k|
			f[4*k] = fo[2*k]; f[4*k + 1] = fo[2*k + 1];
			des[2*k] = deso[k]; dev[2*k] = devo[k];
		};
		(r - 1) do: { |k|
			var aid1, aid2;
			f[4*k + 2] = fo[2*k + 1] + alpha; f[4*k + 3] = fo[2*k + 2] - alpha;
			aid1 = (deso[k] + devo[k]).max(deso[k + 1] + devo[k + 1]);
			aid2 = (deso[k] - devo[k]).min(deso[k + 1] - devo[k + 1]);
			(itype == \type_a).if {
				des[2*k + 1] = (aid1 + aid2)/2;
				dev[2*k + 1] = aid1 - des[2*k + 1];
			} {
				(itype == \type_b).if {
					des[2*k + 1] = 0;
					dev[2*k + 1] = aid1;
				} {
					Error("Type should be either 'type_a' or 'type_b'.").throw
				}
			}
		};
		des = des.stutter(2);
		wt = dev.reciprocal;

		^[f, des, wt]
	}

	/**
	 * Parks-McClellan optimal equiripple FIR order estimator.
	 * [n, fo, ao, w] = FIR.estimateOrder(f, a, dev fs) finds the approximate order n,
	 * normalised frequency band edges fo, frequency band magnitudes ao and weights w
	 * to be used by the FIR.parksMcClellan() function as follows:
	 *    h = FIR.parksMcClellan(n, fo, ao, w)
	 * The resulting filter will approximately meet the specifications given by the
	 * input parameters f, a and dev. f is a vector of cutoff frequencies in Hz,
	 * in ascending order between 0 and half the sampling frequency fs. If you do not
	 * specify fs, it defaults to 2. a is a vector specifying the desired function's
	 * amplitude on the bands defined by f. The length of f is twice the length of a,
	 * minus 2 (it must therefore be even). The first frequency band always starts at
	 * zero, and the last always ends at fs/2. It is not necessary to add these elements
	 * to the f vector. dev is a vector of maximum deviations or ripples (in linear units)
	 * allowable for each band. dev must have the same length as a.
	 *
	 * Example:
	 *    Design a lowpass filter with a passband-edge frequency of 1500Hz, a
	 *    stopband-edge of 2000Hz, passband ripple of 0.01, stopband ripple of 0.1,
	 *    and a sampling frequency of 8000Hz:
	 *
	 *    #n, fo, mo, w = FIR.estimateOrder([1500, 2000], [1, 0], [0.01, 0.1], 8000);
	 *    h = FIR.parksMcClellan(n, fo, mo, w);
	 *
	 * CAUTION 1: The order n is often underestimated. If the filter does not meet the
	 * original specifications, a higher order such as N + 1 or N + 2 will.
	 * CAUTION 2: Results are inaccurate if cutoff frequencies are near zero frequency
	 * or the Nyquist frequency.
	 */
	*estimateOrder { arg fcuts, mags, devs, fsamp = 2;
		var mf, mm, nbands, zz, f1, f2, df, n, remlpord, l, ff, am, aa, wts;

		(fcuts.isNil or: { mags.isNil } or: { devs.isNil }).if {
			Error("Not enough input arguments.").throw
		};
		fcuts = fcuts/fsamp;   // normalise to sampling frequency
		(fcuts any: { |fc| fc > 0.5 }).if {
			Error("Band-edge frequencies must be less than the Nyquist frequency.").throw
		};
		(fcuts any: { |fc| fc < 0 }).if {
			Error("Band-edge frequencies must be positive.").throw
		};

		mf = fcuts.size;
		mm = mags.size;
		nbands = mm;

		(mm != devs.size).if { Error("a and dev must be vectors of the same length.").throw };
		(mf != (2*(nbands - 1))).if { Error("Length of f must be 2*length(a)-2.").throw };

		/**
		 * FIR lowpass filter length estimator
		 * l = remlpord(freq1, freq2, dev1, dev2);
		 *
		 * input:
		 *    freq1: passband cutoff freq (normalised)
		 *    freq2: stopband cutoff freq (normalised)
		 *    dev1:  passband ripple (desired)
		 *    dev2:  stopband attenuation (not in dB)
		 *
		 * output:
		 *    l = filter length (# of samples)
		 *        NOT the order n, which is n = l - 1
		 *
		 * NOTE: Will also work for highpass filters (i.e., f1 > f2)
		 *       Will not work well of transition zone is near f = 0,
		 *       or near f = fs/2
		 */
		remlpord = { |freq1, freq2, delta1, delta2|
			var aa, d1, d2, d, fK, df;
			aa = Matrix.withFlatArray(3, 3, [-4.278e-1, -4.761e-1, 0, -5.941e-1, 7.114e-2, 0, -2.660e-3, 5.309e-3, 0]);
			d1 = log10(delta1);
			d2 = log10(delta2);
			d = Matrix.withFlatArray(1, 3, [1, d1, d1*d1]).mulMatrix(aa).mulMatrix(Matrix.withFlatArray(3, 1, [1, d2, d2*d2]))[0, 0];
			fK = ([1.0, d1 - d2]*[11.01217, 0.51244]).sum;
			df = (freq2 - freq1).abs;
			(d/df - (fK*df) + 1)
		};

		zz = mags collect: { |mag| (mag == 0).asInteger };   // find stopbands
		devs = devs/(zz + mags);                 // divide delta by mag to get relative deviation

		// Determine the smallest width transition zone, separate the passband and stopband edges
		f1 = fcuts[0,2..mf - 2];
		f2 = fcuts[1,3..mf - 1];
		df = (f2 - f1).minItem;
		n = (f2 - f1).minIndex;

		// lowpass case: use formula (ref: Herrmann, Rabiner, Chan)
		(nbands == 2).if {
			l = remlpord.(f1[n], f2[n], devs[0], devs[1]);
		} {
			/**
			 * bandpass case:
			 * - try different lowpasses and take the WORST one that
			 *   goes through the BP specs; try both transition widths
			 * - will also do the bandreject case
			 * - does the multi-band case, one bandpass at a time
			 */
			l = 0;
			(nbands - 2) do: { |i| var l1, l2; i = i + 1;
				l1 = remlpord.(f1[i - 1], f2[i - 1], devs[i], devs[i - 1]);
				l2 = remlpord.(f1[i], f2[i], devs[i], devs[i + 1]);
				l = [l, l1, l2].maxItem;
			}
		};

		n = l.ceil - 1;   // need order, not length for FIR.parksMcClellan()

		ff = [0] ++ (2*fcuts) ++ [1];
		am = 0!(2*nbands - 1);
		(0,2..2*nbands - 2) do: { |i, j| am[i] = mags[j] };
		aa = (am ++ [0]) + ([0] ++ am);
		wts = (1!devs.size)*devs.maxItem/devs;

		// If gain is not zero at Nyquist, the order must be even
		// If the order is odd, we bump up the order by one
		(aa.last != 0 and: { (n%2).booleanValue }).if { n = n + 1 };

		^[n.asInteger, ff, aa, wts]
	}
}