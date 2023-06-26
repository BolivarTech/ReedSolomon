package com.bolivartech.utils.fec.reedsolomon;

/**
 * Copyright 2010,2011,2012,2013 Valhala Networks C.A.<br/>
 *
 * <p>Homepage: <a
 * href="http://www.cuaimacrypt.com">http://www.cuaimacrypt.com</a>.</p>
 * <p>Valhala Networks Homepage: <a
 * href="http://www.valhalanetworks.com">http://www.valhalanetworks.com</a>.</p>
 * <p>
 * This Class is the Valhala Networks's Galois Fields Class.<br/><br/>
 * <p>
 * This class contains utility methods for performing mathematical operations
 * over the Galois Fields. Operations use a given primitive polynomial in
 * calculations.<br/><br/>
 * <p>
 * Throughout this package, elements of the GF are represented as an {@code int}
 * for convenience and speed (but at the cost of memory).<br/><br/>
 * <p>
 * Esta clase implementa utilidades para le manejo de campos de
 * Galois.<br/><br/><br/><br/>
 *
 * @author Julian Bolivar
 * @version 1.0.0
 * @date January 26, 2013.
 * @since 2010, 2011, 2012, 2013
 */
public final class GenericGF {

    /**
     * Campos de Galois estandars
     */
    public static final GenericGF AZTEC_DATA_12 = new GenericGF(0x1069, 4096); // x^12 + x^6 + x^5 + x^3 + 1
    public static final GenericGF AZTEC_DATA_10 = new GenericGF(0x409, 1024); // x^10 + x^3 + 1
    public static final GenericGF AZTEC_DATA_6 = new GenericGF(0x43, 64); // x^6 + x + 1
    public static final GenericGF AZTEC_PARAM = new GenericGF(0x13, 16); // x^4 + x + 1
    public static final GenericGF QR_CODE_FIELD_256 = new GenericGF(0x011D, 256); // x^8 + x^4 + x^3 + x^2 + 1
    public static final GenericGF DATA_MATRIX_FIELD_256 = new GenericGF(0x012D, 256); // x^8 + x^5 + x^3 + x^2 + 1
    public static final GenericGF AZTEC_DATA_8 = DATA_MATRIX_FIELD_256;
    public static final GenericGF MAXICODE_FIELD_64 = AZTEC_DATA_6;

    /*
     * Variables privadas
     */
    private static final int INITIALIZATION_THRESHOLD = 0;
    private final int size;
    private final int primitive;
    private int[] expTable;
    private int[] logTable;
    private GenericGFPoly zero;
    private GenericGFPoly one;
    private boolean initialized = false;

    /**
     * Crea una representacion del campo de Galois de la forma GF(size) usando
     * el polinomio primitivo proporcionado
     *
     * @param primitive Polinomio irreductible cuyos coeficientes son
     *                  representados por los bits de un entero, y en donde los bits memos
     *                  significativos representan el coeficiente constante.
     * @param size      Tamaño de los coeficientes del polinomio
     */
    public GenericGF(int primitive, int size) {
        this.primitive = primitive;
        this.size = size;

        if (size <= INITIALIZATION_THRESHOLD) {
            initialize();
        }
    }

    /**
     * Implementa la suma y la resta (son lo mismo en los campos de Galois(size))
     *
     * @return suma/resta de a y b
     */
    static int addOrSubtract(int a, int b) {
        return a ^ b;
    }

    /*
     * Inicializador del polinomio
     */
    private void initialize() {
        expTable = new int[size];
        logTable = new int[size];
        int x = 1;
        for (int i = 0; i < size; i++) {
            expTable[i] = x;
            x <<= 1; // x = x * 2; we're assuming the generator alpha is 2
            if (x >= size) {
                x ^= primitive;
                x &= size - 1;
            }
        }
        for (int i = 0; i < size - 1; i++) {
            logTable[expTable[i]] = i;
        }
        // logTable[0] == 0 but this should never be used
        zero = new GenericGFPoly(this, new int[]{0});
        one = new GenericGFPoly(this, new int[]{1});
        initialized = true;
    }

    private void checkInit() {
        if (!initialized) {
            initialize();
        }
    }

    GenericGFPoly getZero() {
        checkInit();

        return zero;
    }

    GenericGFPoly getOne() {
        checkInit();

        return one;
    }

    /**
     * @return La representacion monomial de la forma coeficiente * x^grado
     */
    GenericGFPoly buildMonomial(int degree, int coefficient) {
        checkInit();

        if (degree < 0) {
            throw new IllegalArgumentException();
        }
        if (coefficient == 0) {
            return zero;
        }
        int[] coefficients = new int[degree + 1];
        coefficients[0] = coefficient;
        return new GenericGFPoly(this, coefficients);
    }

    /**
     * Implementa la exponenciacion en base 2 en el campo de Galois
     *
     * @return 2 elevado a la potencia de "a" en el GF(size)
     */
    int exp(int a) {
        checkInit();
        return expTable[a];
    }

    /**
     * Implementa el logaritmo en base 2 en el campo de Galois
     *
     * @return Logaritmo en base 2 de "a" en el GF(size)
     */
    int log(int a) {
        checkInit();

        if (a == 0) {
            throw new IllegalArgumentException();
        }
        return logTable[a];
    }

    /**
     * Implementa la inversa multiplicativa en el campo de Galois
     *
     * @return Inversa Multiplicativa de "a"
     */
    int inverse(int a) {
        checkInit();

        if (a == 0) {
            throw new ArithmeticException();
        }
        return expTable[size - logTable[a] - 1];
    }

    /**
     * Implementa la multiplicacion en el campo de Galois
     *
     * @param a Primer coeficiente de la multiplicacion
     * @param b Segundo coeficiente de la multiplicacion
     * @return El producto de "a" y "b" en GF(size)
     */
    int multiply(int a, int b) {
        checkInit();

        if (a == 0 || b == 0) {
            return 0;
        }

        if (a < 0 || b < 0 || a >= size || b >= size) {
            a++;
        }

        int logSum = logTable[a] + logTable[b];
        return expTable[(logSum % size) + logSum / size];
    }

    /**
     * Retorna el tamaño del campo de Galois
     *
     * @return Tamaño del GF
     */
    public int getSize() {
        return size;
    }
}
