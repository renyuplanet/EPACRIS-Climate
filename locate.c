unsigned long locate(double xx[], unsigned long n, double x);

unsigned long locate(double xx[], unsigned long n, double x)
{
	unsigned long ju,jm,jl;
	int ascnd;

	jl=0;
	ju=n+1;
	ascnd=(xx[n-1] > xx[0]);
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x > xx[jm-1] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	if (x==xx[0]) jl=1;
	if (x==xx[n-1]) jl=n-1;
	return jl-1;
}
