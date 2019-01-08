#include <stdio.h>
#include <math.h>
#include <complex.h>
 
double PI;
typedef double complex cplx;

behavior _fft(cplx buf[], cplx out[], int n, int step)
{
	_fft      Feven(out, buf, n, step * 2);
	_fft      Fodd (out + step, buf + step, n, step * 2);
	
	void main(void)
	{
		if (step < n) {
			Feven.main();
			Fodd.main();
	 		int i;
			for (i = 0; i < n; i += 2 * step) {
				cplx t = cexp(-I * PI * i / n) * out[i + step];
				buf[i / 2]     = out[i] + t;
				buf[(i + n)/2] = out[i] - t;
			}
		}
	}
}

behavior fft(cplx buf[], int n){
	_fft     FFT(buf, out, n, 1);
	
	void main(void)
	{
		cplx out[n];
		int i;
		for (i = 0; i < n; i++) out[i] = buf[i];
	 
		FFT.main();
	}
}

behavior print(const char * s, cplx buf[]){
	void main(void) 
	{
		printf("%s", s);
		int i;
		for (i = 0; i < 8; i++)
			if (!cimag(buf[i]))
				printf("%g ", creal(buf[i]));
			else
				printf("(%g, %g) ", creal(buf[i]), cimag(buf[i]));
	}
}
 
behavior Main(void){
	fft    f	   (buf, 8);
	print  printdata("Data: ", buf);
	print  printbuf ("\nFFT : ", buf);
	int main()
	{
		PI = atan2(1, 1) * 4;
		cplx buf[] = {1, 1, 1, 1, 0, 0, 0, 0};
	 	printdata.main();
		f.main();
		printbuf.main();
	 
		return 0;
	}
}
