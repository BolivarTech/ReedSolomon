package com.bolivartech.utils.fec.reedsolomon;

import com.bolivartech.utils.array.ArrayUtils;
import com.bolivartech.utils.converters.Converter;

import java.util.ArrayList;
import java.util.List;

/**
 * Copyright 2010,2011,2012,2013 Valhala Networks C.A.<br/>
 *
 * <p>Homepage: <a
 * href="http://www.cuaimacrypt.com">http://www.cuaimacrypt.com</a>.</p>
 * <p>Valhala Networks Homepage: <a
 * href="http://www.valhalanetworks.com">http://www.valhalanetworks.com</a>.</p>
 * <p>
 * This Class is the Valhala Networks's Reed-Solomon Encoding Class.<br/><br/>
 * <p>
 * Esta clase implementa el metodo Reed-Solomon Encoding.<br/><br/><br/><br/>
 *
 * @author Julian Bolivar
 * @version 1.0.0
 * @date January 26, 2013.
 * @since 2010, 2011, 2012, 2013
 */
public final class ReedSolomonEncoder {

    private final GenericGF field;
    private final List<GenericGFPoly> cachedGenerators;

    /**
     * Constructor por defecto con inicializacion del campo de Galois
     *
     * @param GField Campo de Galois
     */
    public ReedSolomonEncoder(GenericGF GField) {
        if (!GenericGF.QR_CODE_FIELD_256.equals(GField)) {
            throw new IllegalArgumentException("Only QR Code is supported at this time");
        }
        this.field = GField;
        this.cachedGenerators = new ArrayList<GenericGFPoly>();
        cachedGenerators.add(new GenericGFPoly(field, new int[]{1}));
    }

    /**
     * Codifica el mensaje y agrega la cantidad de bytes especificados para la
     * correccion de errores
     *
     * @param Input                Mensaje de entrada original
     * @param ErrorCorrectionBytes Numero de Bytes para la correccion de errores
     * @return Mensaje con los bytes de correccion de errores
     */
    public byte[] encode(byte[] Input, int ErrorCorrectionBytes) {
        int[] toEncode;
        int dataBytes;
        int i;
        int[] infoCoefficients;
        int[] coefficients;
        int numZeroCoefficients;
        byte[] Output;

        if (ErrorCorrectionBytes == 0) {
            throw new IllegalArgumentException("No error correction bytes");
        }
        dataBytes = Input.length;
        if (dataBytes <= 0) {
            throw new IllegalArgumentException("No data bytes provided");
        }
        toEncode = new int[dataBytes + ErrorCorrectionBytes];
        for (i = 0; i < dataBytes; i++) {
            toEncode[i] = Input[i] & 0x000000FF;
        }
        GenericGFPoly generator = buildGenerator(ErrorCorrectionBytes);
        infoCoefficients = new int[dataBytes];
        System.arraycopy(toEncode, 0, infoCoefficients, 0, dataBytes);
        GenericGFPoly info = new GenericGFPoly(field, infoCoefficients);
        info = info.multiplyByMonomial(ErrorCorrectionBytes, 1);
        GenericGFPoly remainder = info.divide(generator)[1];
        coefficients = remainder.getCoefficients();
        numZeroCoefficients = ErrorCorrectionBytes - coefficients.length;
        for (i = 0; i < numZeroCoefficients; i++) {
            toEncode[dataBytes + i] = 0;
        }
        ArrayUtils.arrayCopy(coefficients, 0, toEncode, dataBytes + numZeroCoefficients, coefficients.length);
        Output = Converter.int2byte(toEncode);
        for (i = 0; i < (Output.length / 4); i++) {
            Output[i] = Output[3 + (i * 4)];
        }
        Output = (byte[]) ArrayUtils.resizeArray(Output, Output.length / 4);
        return Output;
    }

    /*
     * Generador de los polinominion del campo de Galois
     */
    private GenericGFPoly buildGenerator(int degree) {
        if (degree >= cachedGenerators.size()) {
            GenericGFPoly lastGenerator = cachedGenerators.get(cachedGenerators.size() - 1);
            for (int d = cachedGenerators.size(); d <= degree; d++) {
                GenericGFPoly nextGenerator = lastGenerator.multiply(new GenericGFPoly(field, new int[]{1, field.exp(d - 1)}));
                cachedGenerators.add(nextGenerator);
                lastGenerator = nextGenerator;
            }
        }
        return cachedGenerators.get(degree);
    }
}
