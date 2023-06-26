package com.valhala.utils.fec.reedsolomon;

/**
 * Copyright 2010,2011,2012,2013 Valhala Networks C.A.<br/>
 *
 * <p>Homepage: <a
 * href="http://www.cuaimacrypt.com">http://www.cuaimacrypt.com</a>.</p>
 * <p>Valhala Networks Homepage: <a
 * href="http://www.valhalanetworks.com">http://www.valhalanetworks.com</a>.</p>
 *
 * This Class is the Valhala Networks's representation of a polynomial whose
 * coefficients are elements of a Galois Fields.<br/> Instances of this class
 * are immutable.<br/><br/>
 *
 * Esta clase implementa la representacion de un polinomio cuyos coeficientes
 * son elemento de un campos de Galois.<br/><br/><br/><br/>
 *
 * @author Julian Bolivar
 * @since 2010,2011,2012,2013
 * @date January 26, 2013.
 * @version 1.0.0
 */
final class GenericGFPoly {

    private final GenericGF field;
    private final int[] coefficients;

    /**
     * Constructor con inicializacion de los coeficientes del polinomio en el
     * campo de Galois.
     *
     * @param GField el la instancia al campo de Galois usado para la efectual
     * los calculos computacionales.
     * @param coefficients Coeficientes del campo de Galois representados como
     * enteros Estos coeficientes estan ordenados desde el mas significativo al
     * menos significativo.
     * @throws IllegalArgumentException si los argumentos son null o vacios o si
     * los coeficientes son 0 o no es una constante polinomial.
     */
    GenericGFPoly(GenericGF GField, int[] coefficients) {
        if (coefficients.length == 0) {
            throw new IllegalArgumentException();
        }
        this.field = GField;
        int coefficientsLength = coefficients.length;
        if (coefficientsLength > 1 && coefficients[0] == 0) {
            // Leading term must be non-zero for anything except the constant polynomial "0"
            int firstNonZero = 1;
            while (firstNonZero < coefficientsLength && coefficients[firstNonZero] == 0) {
                firstNonZero++;
            }
            if (firstNonZero == coefficientsLength) {
                this.coefficients = field.getZero().coefficients;
            } else {
                this.coefficients = new int[coefficientsLength - firstNonZero];
                System.arraycopy(coefficients,
                        firstNonZero,
                        this.coefficients,
                        0,
                        this.coefficients.length);
            }
        } else {
            this.coefficients = coefficients;
        }
    }

    /**
     * Retorna los coeficientes en el polinomio
     *
     * @return Coeficientes del polinomio
     */
    int[] getCoefficients() {
        return coefficients;
    }

    /**
     * Retorna el grado del polinomio
     *
     * @return grado del polinomio
     */
    int getDegree() {
        return coefficients.length - 1;
    }

    /**
     * Retorna verdadero si el polinomio es monomial "0"
     *
     * @return true si el polinomio es monomial "0"
     */
    boolean isZero() {
        return coefficients[0] == 0;
    }

    /**
     * Retorna el coeficiente en el grado x^degree del polinomio.
     *
     * @return Coeficiente del polinomio x^degree
     */
    int getCoefficient(int degree) {
        return coefficients[coefficients.length - 1 - degree];
    }

    /**
     * Retorna el polinomio evaluado en el punto "a"
     *
     * @return El polinomio evaluado en el punto "a"
     */
    int evaluateAt(int a) {
        if (a == 0) {
            // Just return the x^0 coefficient
            return getCoefficient(0);
        }
        int size = coefficients.length;
        if (a == 1) {
            // Just the sum of the coefficients
            int result = 0;
            for (int i = 0; i < size; i++) {
                result = GenericGF.addOrSubtract(result, coefficients[i]);
            }
            return result;
        }
        int result = coefficients[0];
        for (int i = 1; i < size; i++) {
            result = GenericGF.addOrSubtract(field.multiply(a, result), coefficients[i]);
        }
        return result;
    }

    /**
     * Implementa la suma y la resta de los polinomios de Galois
     *
     * @param other El otro polinomio a ser sumado o restado
     * @return Resultado de la operacionde sume o resta
     */
    GenericGFPoly addOrSubtract(GenericGFPoly other) {
        if (!field.equals(other.field)) {
            throw new IllegalArgumentException("GenericGFPolys do not have same GenericGF field");
        }
        if (isZero()) {
            return other;
        }
        if (other.isZero()) {
            return this;
        }

        int[] smallerCoefficients = this.coefficients;
        int[] largerCoefficients = other.coefficients;
        if (smallerCoefficients.length > largerCoefficients.length) {
            int[] temp = smallerCoefficients;
            smallerCoefficients = largerCoefficients;
            largerCoefficients = temp;
        }
        int[] sumDiff = new int[largerCoefficients.length];
        int lengthDiff = largerCoefficients.length - smallerCoefficients.length;
        // Copy high-order terms only found in higher-degree polynomial's coefficients
        System.arraycopy(largerCoefficients, 0, sumDiff, 0, lengthDiff);

        for (int i = lengthDiff; i < largerCoefficients.length; i++) {
            sumDiff[i] = GenericGF.addOrSubtract(smallerCoefficients[i - lengthDiff], largerCoefficients[i]);
        }

        return new GenericGFPoly(field, sumDiff);
    }

    /**
     * Retorna la multilicacion de los polinomios de Galois
     *
     * @param other El otro polinomio a ser multiplicado
     * @return Resultado del polinomio multiplicado
     */
    GenericGFPoly multiply(GenericGFPoly other) {
        if (!field.equals(other.field)) {
            throw new IllegalArgumentException("GenericGFPolys do not have same GenericGF field");
        }
        if (isZero() || other.isZero()) {
            return field.getZero();
        }
        int[] aCoefficients = this.coefficients;
        int aLength = aCoefficients.length;
        int[] bCoefficients = other.coefficients;
        int bLength = bCoefficients.length;
        int[] product = new int[aLength + bLength - 1];
        for (int i = 0; i < aLength; i++) {
            int aCoeff = aCoefficients[i];
            for (int j = 0; j < bLength; j++) {
                product[i + j] = GenericGF.addOrSubtract(product[i + j],
                        field.multiply(aCoeff, bCoefficients[j]));
            }
        }
        return new GenericGFPoly(field, product);
    }

    /**
     * Multiplicacion del polinomio por un escalar
     *
     * @param scalar Valor escalar por el que se va a multiplicar
     * @return Polinomio multiplicado por el escalar
     */
    GenericGFPoly multiply(int scalar) {
        if (scalar == 0) {
            return field.getZero();
        }
        if (scalar == 1) {
            return this;
        }
        int size = coefficients.length;
        int[] product = new int[size];
        for (int i = 0; i < size; i++) {
            product[i] = field.multiply(coefficients[i], scalar);
        }
        return new GenericGFPoly(field, product);
    }

    /**
     * Multiplicacion del polinomio por un monomial
     *
     * @param degree Grado del polinomio
     * @param coefficient Coeficiene monomial
     * @return Multiplicacion por el monomial
     */
    GenericGFPoly multiplyByMonomial(int degree, int coefficient) {
        if (degree < 0) {
            throw new IllegalArgumentException();
        }
        if (coefficient == 0) {
            return field.getZero();
        }
        int size = coefficients.length;
        int[] product = new int[size + degree];
        for (int i = 0; i < size; i++) {
            product[i] = field.multiply(coefficients[i], coefficient);
        }
        return new GenericGFPoly(field, product);
    }

    /**
     * Division del polinomio por otro polinomio
     *
     * @param other Polinomio de Galois por el que vas a dividir
     * @return Polinomio dividido
     */
    GenericGFPoly[] divide(GenericGFPoly other) {
        if (!field.equals(other.field)) {
            throw new IllegalArgumentException("GenericGFPolys do not have same GenericGF field");
        }
        if (other.isZero()) {
            throw new IllegalArgumentException("Divide by 0");
        }

        GenericGFPoly quotient = field.getZero();
        GenericGFPoly remainder = this;

        int denominatorLeadingTerm = other.getCoefficient(other.getDegree());
        int inverseDenominatorLeadingTerm = field.inverse(denominatorLeadingTerm);

        while (remainder.getDegree() >= other.getDegree() && !remainder.isZero()) {
            int degreeDifference = remainder.getDegree() - other.getDegree();
            int scale = field.multiply(remainder.getCoefficient(remainder.getDegree()), inverseDenominatorLeadingTerm);
            GenericGFPoly term = other.multiplyByMonomial(degreeDifference, scale);
            GenericGFPoly iterationQuotient = field.buildMonomial(degreeDifference, scale);
            quotient = quotient.addOrSubtract(iterationQuotient);
            remainder = remainder.addOrSubtract(term);
        }

        return new GenericGFPoly[]{quotient, remainder};
    }

    /**
     * Convierte el polinomio en su representacion de String
     *
     * @return Representacion en String del Polinomio
     */
    @Override
    public String toString() {
        StringBuilder result = new StringBuilder(8 * getDegree());
        for (int degree = getDegree(); degree >= 0; degree--) {
            int coefficient = getCoefficient(degree);
            if (coefficient != 0) {
                if (coefficient < 0) {
                    result.append(" - ");
                    coefficient = -coefficient;
                } else {
                    if (result.length() > 0) {
                        result.append(" + ");
                    }
                }
                if (degree == 0 || coefficient != 1) {
                    int alphaPower = field.log(coefficient);
                    if (alphaPower == 0) {
                        result.append('1');
                    } else if (alphaPower == 1) {
                        result.append('a');
                    } else {
                        result.append("a^");
                        result.append(alphaPower);
                    }
                }
                if (degree != 0) {
                    if (degree == 1) {
                        result.append('x');
                    } else {
                        result.append("x^");
                        result.append(degree);
                    }
                }
            }
        }
        return result.toString();
    }
}
