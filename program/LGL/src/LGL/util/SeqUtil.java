/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package LGL.util;

/**
 *
 * @author ligl
 */
public class SeqUtil {

    final static char FORWARD = '+';
    final static char REVERSE = '-';
    private static char[] complementTable = new char[255];

    static {
        complementTable['A'] = 'T';
        complementTable['C'] = 'G';
        complementTable['G'] = 'C';
        complementTable['T'] = 'A';
        complementTable['N'] = 'N';

        complementTable['a'] = 't';
        complementTable['c'] = 'g';
        complementTable['g'] = 'c';
        complementTable['t'] = 'a';
        complementTable['n'] = 'n';
    }

    public static boolean isForwardStrand(char strand) {
        if ((strand == 'F') || (strand == 'f') || (strand == '+')) {
            return true;
        } else {
            return false;
        }
    }

    public static int strand2num(char strand) {
        if (isForwardStrand(strand)) {
            return 1;
        } else {
            return 0;
        }
    }

    public static char num2strand(int strandId) {
        if (strandId == 1) {
            return FORWARD;
        } else {
            return REVERSE;
        }
    }

    public static char getStrand(char strand) {
        if ((strand == 'F') || (strand == 'f') || (strand == '+')) {
            return '+';
        } else if ((strand == 'R') || (strand == 'r') || (strand == '-')) {
            return '-';
        } else {
            System.out.println("The strand is invalid: " + strand);
            System.out.println("Replaced the strand with \".\"");
            return '.';
        }
    }

    public static char InvertStrand(char strand) {
        if ((strand == 'F') || (strand == 'f') || (strand == '+')) {
            return '-';
        } else if ((strand == 'R') || (strand == 'r') || (strand == '-')) {
            return '+';
        } else {
            System.out.println("The strand is invalid: " + strand);
            System.out.println("Replaced the strand with \".\"");
            return '.';
        }
    }

    public static String revComplement(String seq) {
        StringBuilder result = new StringBuilder(seq);
        result.reverse();
        for (int i = seq.length() - 1; i >= 0; i--) {
            switch (result.charAt(i)) {
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'N':
                case 'a':
                case 'c':
                case 'g':
                case 't':
                case 'n':
                    result.setCharAt(i, complementTable[result.charAt(i)]);
                    break;
                default:
                    break;
            }
        }
        return result.toString();
    }

    public static String reverse(String seq) {
        StringBuilder result = new StringBuilder(seq);
        result.reverse();
        return result.toString();
    }

    public static String complement(String seq) {
        StringBuilder result = new StringBuilder(seq);
        for (int i = seq.length() - 1; i >= 0; i--) {
            switch (result.charAt(i)) {
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'N':
                case 'a':
                case 'c':
                case 'g':
                case 't':
                case 'n':
                    result.setCharAt(i, complementTable[result.charAt(i)]);
                    break;
                default:
                    break;
            }
        }
        return result.toString();
    }

    public static void main(String[] args) {
        System.out.println("AGGTCCAGTATGATTTGATTGGTTGGATCATATATC");
        System.out.println(revComplement("AGGTCCAGTATGATTTGATTGGTTGGATCATATATC"));
        System.out.println(revComplement(revComplement("AGGTCCAGTATGATTTGATTGGTTGGATCATATATC")));

        if (args.length >= 1) {
            System.out.println("");
            System.out.println("original:           " + args[0]);
            System.out.println("reverse:            " + reverse(args[0]));
            System.out.println("complement:         " + complement(args[0]));
            System.out.println("reverse complement: " + revComplement(args[0]));
        }
    }
}
