template <class T>
struct F1dim {
	const VecDoub &p;
	const VecDoub &xi;
	Int n;
	T &func;
	VecDoub xt;
	F1dim(VecDoub_I &pp, VecDoub_I &xii, T &funcc) : p(pp),
		xi(xii), n(pp.size()), func(funcc), xt(n) {}
	Doub operator() (const Doub x)
	{
		for (Int j = 0;j<n;j++)
			xt[j] = p[j] + x*xi[j];
		return func(xt);
	}
};
template <class T>
struct Linemethod {
	VecDoub p;
	VecDoub xi;
	T &func;
	Int n;
	Linemethod(T &funcc) : func(funcc) {}
	Doub linmin()
	{
		Doub ax, xx, xmin;
		n = p.size();
		F1dim<T> f1dim(p, xi, func);
		ax = 0.0;
		xx = 1.0;
		Brent brent;
		brent.bracket(ax, xx, f1dim);
		xmin = brent.minimize(f1dim);
		for (Int j = 0;j<n;j++) {
			xi[j] *= xmin;
			p[j] += xi[j];
		}
		return brent.fmin;
	}
};
template <class T>
struct Df1dim {
	const VecDoub &p;
	const VecDoub &xi;
	Int n;
	T &funcd;
	VecDoub xt;
	VecDoub dft;
	Df1dim(VecDoub_I &pp, VecDoub_I &xii, T &funcdd) : p(pp),
		xi(xii), n(pp.size()), funcd(funcdd), xt(n), dft(n) {}
	Doub operator()(const Doub x)
	{
		for (Int j = 0;j<n;j++)
			xt[j] = p[j] + x*xi[j];
		return funcd(xt);
	}
	Doub df(const Doub x)
	{
		Doub df1 = 0.0;
		funcd.df(xt, dft);
		for (Int j = 0;j<n;j++)
			df1 += dft[j] * xi[j];
		return df1;
	}
};
template <class T>
struct Dlinemethod {
	VecDoub p;
	VecDoub xi;
	T &func;
	Int n;
	Dlinemethod(T &funcc) : func(funcc) {}
	Doub linmin()
	{
		Doub ax, xx, xmin;
		n = p.size();
		Df1dim<T> df1dim(p, xi, func);
		ax = 0.0;
		xx = 1.0;
		Dbrent dbrent;
		dbrent.bracket(ax, xx, df1dim);
		xmin = dbrent.minimize(df1dim);
		for (Int j = 0;j<n;j++) {
			xi[j] *= xmin;
			p[j] += xi[j];
		}
		return dbrent.fmin;
	}
};
template <class T>
struct Powell : Linemethod<T> {
	Int iter;
	Doub fret;
	using Linemethod<T>::func;
	using Linemethod<T>::linmin;
	using Linemethod<T>::p;
	using Linemethod<T>::xi;
	const Doub ftol;
	Powell(T &func, const Doub ftoll = 3.0e-8) : Linemethod<T>(func),
		ftol(ftoll) {}
	VecDoub minimize(VecDoub_I &pp)
	{
		Int n = pp.size();
		MatDoub ximat(n, n, 0.0);
		for (Int i = 0;i<n;i++) ximat[i][i] = 1.0;
		return minimize(pp, ximat);
	}
	VecDoub minimize(VecDoub_I &pp, MatDoub_IO &ximat)
	{
		const Int ITMAX = 200;
		const Doub TINY = 1.0e-25;
		Doub fptt;
		Int n = pp.size();
		p = pp;
		VecDoub pt(n), ptt(n);
		xi.resize(n);
		fret = func(p);
		for (Int j = 0;j<n;j++) pt[j] = p[j];
		for (iter = 0;;++iter) {
			Doub fp = fret;
			Int ibig = 0;
			Doub del = 0.0;
			for (Int i = 0;i<n;i++) {
				for (Int j = 0;j<n;j++) xi[j] = ximat[j][i];
				fptt = fret;
				fret = linmin();
				if (fptt - fret > del) {
					del = fptt - fret;
					ibig = i + 1;
				}
			}
			if (2.0*(fp - fret) <= ftol*(abs(fp) + abs(fret)) + TINY) {
				return p;
			}
			//if (iter == ITMAX) throw("powell exceeding maximum iterations.");
			if (iter == ITMAX) {
				break;
			}
			for (Int j = 0;j<n;j++) {
				ptt[j] = 2.0*p[j] - pt[j];
				xi[j] = p[j] - pt[j];
				pt[j] = p[j];
			}
			fptt = func(ptt);
			if (fptt < fp) {
				Doub t = 2.0*(fp - 2.0*fret + fptt)*SQR(fp - fret - del) - del*SQR(fp - fptt);
				if (t < 0.0) {
					fret = linmin();
					for (Int j = 0;j<n;j++) {
						ximat[j][ibig - 1] = ximat[j][n - 1];
						ximat[j][n - 1] = xi[j];
					}
				}
			}
		}
	}
};
template <class T>
struct Frprmn : Linemethod<T> {
	Int iter;
	Doub fret;
	using Linemethod<T>::func;
	using Linemethod<T>::linmin;
	using Linemethod<T>::p;
	using Linemethod<T>::xi;
	const Doub ftol;
	Frprmn(T &funcd, const Doub ftoll = 3.0e-8) : Linemethod<T>(funcd),
		ftol(ftoll) {}
	VecDoub minimize(VecDoub_I &pp)
	{
		const Int ITMAX = 200;
		const Doub EPS = 1.0e-18;
		const Doub GTOL = 1.0e-8;
		Doub gg, dgg;
		Int n = pp.size();
		p = pp;
		VecDoub g(n), h(n);
		xi.resize(n);
		Doub fp = func(p);
		func.df(p, xi);
		for (Int j = 0;j<n;j++) {
			g[j] = -xi[j];
			xi[j] = h[j] = g[j];
		}
		for (Int its = 0;its<ITMAX;its++) {
			iter = its;
			fret = linmin();
			if (2.0*abs(fret - fp) <= ftol*(abs(fret) + abs(fp) + EPS))
				return p;
			fp = fret;
			func.df(p, xi);
			Doub test = 0.0;
			Doub den = MAX(abs(fp), 1.0);
			for (Int j = 0;j<n;j++) {
				Doub temp = abs(xi[j])*MAX(abs(p[j]), 1.0) / den;
				if (temp > test) test = temp;
			}
			if (test < GTOL) return p;
			dgg = gg = 0.0;
			for (Int j = 0;j<n;j++) {
				gg += g[j] * g[j];
				//			  dgg += xi[j]*xi[j];
				dgg += (xi[j] + g[j])*xi[j];
			}
			if (gg == 0.0)
				return p;
			Doub gam = dgg / gg;
			for (Int j = 0;j<n;j++) {
				g[j] = -xi[j];
				xi[j] = h[j] = g[j] + gam*h[j];
			}
		}
		throw("Too many iterations in frprmn");
	}
};

struct Bracketmethod {
	Doub ax, bx, cx, fa, fb, fc;
	template <class T>
	void bracket(const Doub a, const Doub b, T &func)
	{
		const Doub GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
		ax = a; bx = b;
		Doub fu;
		fa = func(ax);
		fb = func(bx);
		if (fb > fa) {
			SWAP(ax, bx);
			SWAP(fb, fa);
		}
		cx = bx + GOLD*(bx - ax);
		fc = func(cx);
		while (fb > fc) {
			Doub r = (bx - ax)*(fb - fc);
			Doub q = (bx - cx)*(fb - fa);
			Doub u = bx - ((bx - cx)*q - (bx - ax)*r) /
				(2.0*SIGN(MAX(abs(q - r), TINY), q - r));
			Doub ulim = bx + GLIMIT*(cx - bx);
			if ((bx - u)*(u - cx) > 0.0) {
				fu = func(u);
				if (fu < fc) {
					ax = bx;
					bx = u;
					fa = fb;
					fb = fu;
					return;
				}
				else if (fu > fb) {
					cx = u;
					fc = fu;
					return;
				}
				u = cx + GOLD*(cx - bx);
				fu = func(u);
			}
			else if ((cx - u)*(u - ulim) > 0.0) {
				fu = func(u);
				if (fu < fc) {
					shft3(bx, cx, u, u + GOLD*(u - cx));
					shft3(fb, fc, fu, func(u));
				}
			}
			else if ((u - ulim)*(ulim - cx) >= 0.0) {
				u = ulim;
				fu = func(u);
			}
			else {
				u = cx + GOLD*(cx - bx);
				fu = func(u);
			}
			shft3(ax, bx, cx, u);
			shft3(fa, fb, fc, fu);
		}
	}
	inline void shft2(Doub &a, Doub &b, const Doub c)
	{
		a = b;
		b = c;
	}
	inline void shft3(Doub &a, Doub &b, Doub &c, const Doub d)
	{
		a = b;
		b = c;
		c = d;
	}
	inline void mov3(Doub &a, Doub &b, Doub &c, const Doub d, const Doub e,
		const Doub f)
	{
		a = d; b = e; c = f;
	}
};
struct Golden : Bracketmethod {
	Doub xmin, fmin;
	const Doub tol;
	Golden(const Doub toll = 3.0e-8) : tol(toll) {}
	template <class T>
	Doub minimize(T &func)
	{
		const Doub R = 0.61803399, C = 1.0 - R;
		Doub x1, x2;
		Doub x0 = ax;
		Doub x3 = cx;
		if (abs(cx - bx) > abs(bx - ax)) {
			x1 = bx;
			x2 = bx + C*(cx - bx);
		}
		else {
			x2 = bx;
			x1 = bx - C*(bx - ax);
		}
		Doub f1 = func(x1);
		Doub f2 = func(x2);
		while (abs(x3 - x0) > tol*(abs(x1) + abs(x2))) {
			if (f2 < f1) {
				shft3(x0, x1, x2, R*x2 + C*x3);
				shft2(f1, f2, func(x2));
			}
			else {
				shft3(x3, x2, x1, R*x1 + C*x0);
				shft2(f2, f1, func(x1));
			}
		}
		if (f1 < f2) {
			xmin = x1;
			fmin = f1;
		}
		else {
			xmin = x2;
			fmin = f2;
		}
		return xmin;
	}
};
struct Brent : Bracketmethod {
	Doub xmin, fmin;
	const Doub tol;
	Brent(const Doub toll = 3.0e-8) : tol(toll) {}
	template <class T>
	Doub minimize(T &func)
	{
		const Int ITMAX = 100;
		const Doub CGOLD = 0.3819660;
		const Doub ZEPS = numeric_limits<Doub>::epsilon()*1.0e-3;
		Doub a, b, d = 0.0, etemp, fu, fv, fw, fx;
		Doub p, q, r, tol1, tol2, u, v, w, x, xm;
		Doub e = 0.0;

		a = (ax < cx ? ax : cx);
		b = (ax > cx ? ax : cx);
		x = w = v = bx;
		fw = fv = fx = func(x);
		for (Int iter = 0;iter<ITMAX;iter++) {
			xm = 0.5*(a + b);
			tol2 = 2.0*(tol1 = tol*abs(x) + ZEPS);
			if (abs(x - xm) <= (tol2 - 0.5*(b - a))) {
				fmin = fx;
				return xmin = x;
			}
			if (abs(e) > tol1) {
				r = (x - w)*(fx - fv);
				q = (x - v)*(fx - fw);
				p = (x - v)*q - (x - w)*r;
				q = 2.0*(q - r);
				if (q > 0.0) p = -p;
				q = abs(q);
				etemp = e;
				e = d;
				if (abs(p) >= abs(0.5*q*etemp) || p <= q*(a - x)
					|| p >= q*(b - x))
					d = CGOLD*(e = (x >= xm ? a - x : b - x));
				else {
					d = p / q;
					u = x + d;
					if (u - a < tol2 || b - u < tol2)
						d = SIGN(tol1, xm - x);
				}
			}
			else {
				d = CGOLD*(e = (x >= xm ? a - x : b - x));
			}
			u = (abs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
			fu = func(u);
			if (fu <= fx) {
				if (u >= x) a = x; else b = x;
				shft3(v, w, x, u);
				shft3(fv, fw, fx, fu);
			}
			else {
				if (u < x) a = u; else b = u;
				if (fu <= fw || w == x) {
					v = w;
					w = u;
					fv = fw;
					fw = fu;
				}
				else if (fu <= fv || v == x || v == w) {
					v = u;
					fv = fu;
				}
			}
		}
		throw("Too many iterations in brent");
	}
};
struct Dbrent : Bracketmethod {
	Doub xmin, fmin;
	const Doub tol;
	Dbrent(const Doub toll = 3.0e-8) : tol(toll) {}
	template <class T>
	Doub minimize(T &funcd)
	{
		const Int ITMAX = 100;
		const Doub ZEPS = numeric_limits<Doub>::epsilon()*1.0e-3;
		Bool ok1, ok2;
		Doub a, b, d = 0.0, d1, d2, du, dv, dw, dx, e = 0.0;
		Doub fu, fv, fw, fx, olde, tol1, tol2, u, u1, u2, v, w, x, xm;

		a = (ax < cx ? ax : cx);
		b = (ax > cx ? ax : cx);
		x = w = v = bx;
		fw = fv = fx = funcd(x);
		dw = dv = dx = funcd.df(x);
		for (Int iter = 0;iter<ITMAX;iter++) {
			xm = 0.5*(a + b);
			tol1 = tol*abs(x) + ZEPS;
			tol2 = 2.0*tol1;
			if (abs(x - xm) <= (tol2 - 0.5*(b - a))) {
				fmin = fx;
				return xmin = x;
			}
			if (abs(e) > tol1) {
				d1 = 2.0*(b - a);
				d2 = d1;
				if (dw != dx) d1 = (w - x)*dx / (dx - dw);
				if (dv != dx) d2 = (v - x)*dx / (dx - dv);
				u1 = x + d1;
				u2 = x + d2;
				ok1 = (a - u1)*(u1 - b) > 0.0 && dx*d1 <= 0.0;
				ok2 = (a - u2)*(u2 - b) > 0.0 && dx*d2 <= 0.0;
				olde = e;
				e = d;
				if (ok1 || ok2) {
					if (ok1 && ok2)
						d = (abs(d1) < abs(d2) ? d1 : d2);
					else if (ok1)
						d = d1;
					else
						d = d2;
					if (abs(d) <= abs(0.5*olde)) {
						u = x + d;
						if (u - a < tol2 || b - u < tol2)
							d = SIGN(tol1, xm - x);
					}
					else {
						d = 0.5*(e = (dx >= 0.0 ? a - x : b - x));
					}
				}
				else {
					d = 0.5*(e = (dx >= 0.0 ? a - x : b - x));
				}
			}
			else {
				d = 0.5*(e = (dx >= 0.0 ? a - x : b - x));
			}
			if (abs(d) >= tol1) {
				u = x + d;
				fu = funcd(u);
			}
			else {
				u = x + SIGN(tol1, d);
				fu = funcd(u);
				if (fu > fx) {
					fmin = fx;
					return xmin = x;
				}
			}
			du = funcd.df(u);
			if (fu <= fx) {
				if (u >= x) a = x; else b = x;
				mov3(v, fv, dv, w, fw, dw);
				mov3(w, fw, dw, x, fx, dx);
				mov3(x, fx, dx, u, fu, du);
			}
			else {
				if (u < x) a = u; else b = u;
				if (fu <= fw || w == x) {
					mov3(v, fv, dv, w, fw, dw);
					mov3(w, fw, dw, u, fu, du);
				}
				else if (fu < fv || v == x || v == w) {
					mov3(v, fv, dv, u, fu, du);
				}
			}
		}
		throw("Too many iterations in routine dbrent");
	}
};
