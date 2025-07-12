

package com.harium.keel.catalano.math;


public class ComplexNumber {

  
    public double real = 0;
  
    public double imaginary = 0;

  
    public ComplexNumber() {
        this(0, 0);
    }

   
    public ComplexNumber(double real, double imaginary) {
        this.real = real;
        this.imaginary = imaginary;
    }

   
    public ComplexNumber(ComplexNumber z1) {
        this.real = z1.real;
        this.imaginary = z1.imaginary;
    }

    public double getMagnitude() {
        return Math.sqrt(real * real + imaginary * imaginary);
    }

   
    public double getSquaredMagnitude() {
        return real * real + imaginary * imaginary;
    }

    
    public double getPhase() {
        return Math.atan2(imaginary, real);
    }

  
    public static double[] getReal(ComplexNumber[] cn) {
        double[] n = new double[cn.length];
        for (int i = 0; i < n.length; i++) {
            n[i] = cn[i].real;
        }
        return n;
    }


    public static double[] getImaginary(ComplexNumber[] cn) {
        double[] n = new double[cn.length];
        for (int i = 0; i < n.length; i++) {
            n[i] = cn[i].imaginary;
        }
        return n;
    }

   
    public static double[][] getReal(ComplexNumber[][] cn) {
        double[][] n = new double[cn.length][cn[0].length];
        for (int i = 0; i < n.length; i++) {
            for (int j = 0; j < n[0].length; j++) {
                n[i][j] = cn[i][j].real;
            }
        }
        return n;
    }

  
    public static double[][] getImaginary(ComplexNumber[][] cn) {
        double[][] n = new double[cn.length][cn[0].length];
        for (int i = 0; i < n.length; i++) {
            for (int j = 0; j < n[0].length; j++) {
                n[i][j] = cn[i][j].imaginary;
            }
        }
        return n;
    }

    
    public static void Swap(ComplexNumber z1) {
        double t = z1.real;
        z1.real = z1.imaginary;
        z1.imaginary = t;
    }

  
    public static void Swap(ComplexNumber[] z) {
        for (int i = 0; i < z.length; i++) {
            z[i] = new ComplexNumber(z[i].imaginary, z[i].real);
        }
    }

   
    public static void Swap(ComplexNumber[][] z) {
        for (int i = 0; i < z.length; i++) {
            for (int j = 0; j < z[0].length; j++) {
                z[i][j] = new ComplexNumber(z[i][j].imaginary, z[i][j].real);
            }
        }
    }

  
    public static double Abs(ComplexNumber z) {
        return Magnitude(z);
    }

  
    public static double[] Abs(ComplexNumber[] z) {
        double[] values = new double[z.length];
        for (int i = 0; i < values.length; i++) {
            values[i] = z[i].getMagnitude();
        }
        return values;
    }

  
    public static double[][] Abs(ComplexNumber[][] z) {
        double[][] values = new double[z.length][z[0].length];
        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < values[0].length; j++) {
                values[i][j] = z[i][j].getMagnitude();
            }
        }
        return values;
    }

  
    public static ComplexNumber Add(ComplexNumber z1, ComplexNumber z2) {
        return new ComplexNumber(z1.real + z2.real, z1.imaginary + z2.imaginary);
    }

  
    public static ComplexNumber Add(ComplexNumber z1, double scalar) {
        return new ComplexNumber(z1.real + scalar, z1.imaginary);
    }

 
    public void Add(double scalar) {
        this.real += scalar;
    }

 
    public static ComplexNumber Subtract(ComplexNumber z1, ComplexNumber z2) {
        return new ComplexNumber(z1.real - z2.real, z1.imaginary - z2.imaginary);
    }

 
    public static ComplexNumber Subtract(ComplexNumber z1, double scalar) {
        return new ComplexNumber(z1.real - scalar, z1.imaginary);
    }

    
    public void Subtract(double scalar) {
        this.real -= scalar;
    }

  
    public static double Magnitude(ComplexNumber z) {
        return Math.sqrt(z.real * z.real + z.imaginary * z.imaginary);
    }


    public static ComplexNumber Multiply(ComplexNumber z1, ComplexNumber z2) {
        double z1R = z1.real, z1I = z1.imaginary;
        double z2R = z2.real, z2I = z2.imaginary;

        return new ComplexNumber(z1R * z2R - z1I * z2I, z1R * z2I + z1I * z2R);
    }

  
    public static ComplexNumber Multiply(ComplexNumber z1, double scalar) {
        return new ComplexNumber(z1.real * scalar, z1.imaginary * scalar);
    }

  
    public void Multiply(double scalar) {
        this.real *= scalar;
        this.imaginary *= scalar;
    }

 
    public static ComplexNumber Divide(ComplexNumber z1, ComplexNumber z2) {

        ComplexNumber conj = ComplexNumber.Conjugate(z2);

        double a = z1.real * conj.real + ((z1.imaginary * conj.imaginary) * -1);
        double b = z1.real * conj.imaginary + (z1.imaginary * conj.real);

        double c = z2.real * conj.real + ((z2.imaginary * conj.imaginary) * -1);

        return new ComplexNumber(a / c, b / c);
    }

  
    public void Divide(ComplexNumber z1) {
        ComplexNumber conj = ComplexNumber.Conjugate(z1);

        double a = this.real * conj.real + ((this.imaginary * conj.imaginary) * -1);
        double b = this.real * conj.imaginary + (this.imaginary * conj.real);

        double c = z1.real * conj.real + ((z1.imaginary * conj.imaginary) * -1);

        this.real = a / c;
        this.imaginary = b / c;
    }

  
    public static ComplexNumber Divide(ComplexNumber z1, double scalar) {
        return new ComplexNumber(z1.real / scalar, z1.imaginary / scalar);
    }

  
    public void Divide(double scalar) {

        if (scalar == 0) {
            try {
                throw new ArithmeticException("Can not divide by zero.");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        this.real /= scalar;
        this.imaginary /= scalar;
    }


    public static ComplexNumber Pow(ComplexNumber z1, double n) {

        double norm = Math.pow(z1.getMagnitude(), n);
        double angle = 360 - Math.abs(Math.toDegrees(Math.atan(z1.imaginary / z1.real)));

        double common = n * angle;

        double r = norm * Math.cos(Math.toRadians(common));
        double i = norm * Math.sin(Math.toRadians(common));

        return new ComplexNumber(r, i);

    }

  
    public void Pow(double n) {
        double norm = Math.pow(getMagnitude(), n);
        double angle = 360 - Math.abs(Math.toDegrees(Math.atan(this.imaginary / this.real)));

        double common = n * angle;

        this.real = norm * Math.cos(Math.toRadians(common));
        this.imaginary = norm * Math.sin(Math.toRadians(common));
    }

  
    public static ComplexNumber Log(ComplexNumber z1) {
        ComplexNumber result = new ComplexNumber();

        if ((z1.real > 0.0) && (z1.imaginary == 0.0)) {
            result.real = Math.log(z1.real);
            result.imaginary = 0.0;
        } else if (z1.real == 0.0) {
            if (z1.imaginary > 0.0) {
                result.real = Math.log(z1.imaginary);
                result.imaginary = Math.PI / 2.0;
            } else {
                result.real = Math.log(-(z1.imaginary));
                result.imaginary = -Math.PI / 2.0;
            }
        } else {
            result.real = Math.log(z1.getMagnitude());
            result.imaginary = Math.atan2(z1.imaginary, z1.real);
        }

        return result;
    }

 
    public static ComplexNumber Exp(ComplexNumber z1) {
        ComplexNumber x, y;
        x = new ComplexNumber(Math.exp(z1.real), 0.0);
        y = new ComplexNumber(Math.cos(z1.imaginary), Math.sin(z1.imaginary));

        return Multiply(x, y);
    }

 
    public static ComplexNumber Sin(ComplexNumber z1) {
        ComplexNumber result = new ComplexNumber();

        if (z1.imaginary == 0.0) {
            result.real = Math.sin(z1.real);
            result.imaginary = 0.0;
        } else {
            result.real = Math.sin(z1.real) * Math.cosh(z1.imaginary);
            result.imaginary = Math.cos(z1.real) * Math.sinh(z1.imaginary);
        }

        return result;
    }


    public static ComplexNumber Cos(ComplexNumber z1) {
        ComplexNumber result = new ComplexNumber();

        if (z1.imaginary == 0.0) {
            result.real = Math.cos(z1.real);
            result.imaginary = 0.0;
        } else {
            result.real = Math.cos(z1.real) * Math.cosh(z1.imaginary);
            result.imaginary = -Math.sin(z1.real) * Math.sinh(z1.imaginary);
        }

        return result;
    }

   
    public static ComplexNumber Tan(ComplexNumber z1) {
        ComplexNumber result = new ComplexNumber();

        if (z1.imaginary == 0.0) {
            result.real = Math.tan(z1.real);
            result.imaginary = 0.0;
        } else {
            double real2 = 2 * z1.real;
            double imag2 = 2 * z1.imaginary;
            double denom = Math.cos(real2) + Math.cosh(real2);

            result.real = Math.sin(real2) / denom;
            result.imaginary = Math.sinh(imag2) / denom;
        }

        return result;
    }

   
    public void Conjugate() {
        this.imaginary *= -1;
    }


    public static ComplexNumber Conjugate(ComplexNumber z1) {
        return new ComplexNumber(z1.real, z1.imaginary * -1);
    }

    @Override
    public String toString() {
        if (this.imaginary >= 0)
            return this.real + " +" + this.imaginary + "i";
        return this.real + " " + this.imaginary + "i";
    }
}
