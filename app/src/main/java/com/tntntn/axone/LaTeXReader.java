package com.tntntn.axone;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

public class LaTeXReader {
    public static List<Lexeme> latexPrepare(String input) {
        input = input.replace("\\left(", "(")
                .replace("\\right)", ")")
                .replace("Gamma^{-1}", "invGamma")
                .replace("sin^{-1}", "arcsin")
                .replace("cos^{-1}", "arccos")
                .replace("tan^{-1}", "arctan")
                .replace("cot^{-1}", "arccot")
                .replace("sec^{-1}", "arcsec")
                .replace("csc^{-1}", "arccsc")
                .replace("tg^{-1}", "arctan")
                .replace("ctg^{-1}", "arccot")
                .replace("cosec^{-1}", "arccsc")
                .replace("sinh^{-1}", "arcsinh")
                .replace("sh^{-1}", "arcsinh")
                .replace("cosh^{-1}", "arccosh")
                .replace("ch^{-1}", "arccosh")
                .replace("tanh^{-1}", "arctanh")
                .replace("th^{-1}", "arctanh")
                .replace("coth^{-1}", "arccoth")
                .replace("cth^{-1}", "arccoth")
                .replace("sech^{-1}", "arcsech")
                .replace("sch^{-1}", "arcsech")
                .replace("csch^{-1}", "arccsch")
                .replace("sinp^{-1}", "arcsinp")
                .replace("cosp^{-1}", "arccosp")
                .replace("tanp^{-1}", "arctanp")
                .replace("cotp^{-1}", "arccotp")
                .replace("secp^{-1}", "arcsecp")
                .replace("cscp^{-1}", "arccscp")
                .replace(")(", ")*(")
                .replace("1_p", "isPrime")
                .replace("1_P", "isPrime")
                .replace("1_{p}", "isPrime")
                .replace("1_{P}", "isPrime")
                .replace("1_{\\mathbb{p}}", "isPrime")
                .replace("1_{\\mathbb{P}}", "isPrime")
                .replace("exp^infty", "expinf")
                .replace("exp^inf", "expinf")
                .replace("exp^{inf}", "expinf")
                .replace("exp^{infty}", "expinf");
        int i = 0;
        int to = input.length();
        while (i < to) {
            if (i != to - 1 && input.charAt(i) == ')' && (i<to-1?(Character.isAlphabetic(input.charAt(i + 1)) || Character.isDigit(input.charAt(i + 1))):false)) {
                input = input.substring(0, i+1) + "*" + input.substring(i+1);
                to++;
            }
            if (Character.isDigit(input.charAt(i)) && (i<to-1?Character.isAlphabetic(input.charAt(i + 1)):false)) {
                input = input.substring(0, i+1) + "*" + input.substring(i+1);
                to++;
            }
            if (input.charAt(i) == '(' && (i > 0 ? Character.isDigit(input.charAt(i - 1)) : false) && (i > 1 ? input.charAt(i - 2) != '_' && input.charAt(i - 2) != '{' : true)) {
                input = input.substring(0, i) + "*" + input.substring(i);
                to++;
            }
            //if (input.charAt(i) == '(' && (i > 0 ? Character.isAlphabetic(input.charAt(i - 1)) : false)) {
            //    input = input.substring(0, i) + "*" + input.substring(i);
            //    to++;
            //}
            i++;
        }
        input = MainDraw.pow_converter(input, 0);
        List<Lexeme> lexemes = lexAnalyze(input);
        i = 0;
        to = lexemes.size()-1;
        while (i < to) {
            if (lexemes.get(i).type==LexemeType.RIGHT_BRACKET && lexemes.get(i+1).type==LexemeType.LEFT_BRACKET) {
                to++;
                lexemes.add(i+1, new Lexeme(LexemeType.OP_MUL,"*"));
            }
            i++;
        }
        return lexemes;
    }

    public static double latexCalc(List<Lexeme> lexemes, double x, double y) {
        functionMap.put("x", args -> x * (args.isEmpty() ? 1 : args.get(0)));
        functionMap.put("y", args -> y * (args.isEmpty() ? 1 : args.get(0)));
        LexemeBuffer lexemeBuffer = new LexemeBuffer(lexemes);
        return expr(lexemeBuffer);
    }

    public enum LexemeType {
        LEFT_BRACKET, RIGHT_BRACKET,
        OP_PLUS, OP_MINUS, OP_MUL, OP_DIV,
        NUMBER, NAME, COMMA,
        EOF
    }

    public interface Function {
        double apply(List<Double> args);
    }

    private static final HashMap<String, Function> functionMap;
    private static final Set<String> constFunctions;

    static {
        HashMap<String, Function> constFunctionMap = new HashMap<String, Function>() {{
            put("qpi", args -> Math.PI / 4 * (args.isEmpty() ? 1 : args.get(0)));
            put("hpi", args -> Math.PI / 2 * (args.isEmpty() ? 1 : args.get(0)));
            put("tau", args -> Math.PI * 2 * (args.isEmpty() ? 1 : args.get(0)));
            put("pi", args -> Math.PI * (args.isEmpty() ? 1 : args.get(0)));
            put("e", args -> Math.E * (args.isEmpty() ? 1 : args.get(0)));
            put("phi", args -> Functions.GOLDEN_RATIO * (args.isEmpty() ? 1 : args.get(0)));
            put("D", args -> Functions.DOTTIE * (args.isEmpty() ? 1 : args.get(0)));
            put("x", null);
            put("y", null);
        }};
        constFunctions = constFunctionMap.keySet();

        functionMap = new HashMap<String, Function>() {{
            put("sgn", args -> {if (args.size()==1) {return Functions.sgn(args.get(0));} return Double.NaN;});
            put("mod", args -> {if (args.size()==2) {return Functions.mod(args.get(0),args.get(1));} return Double.NaN;});
            put("sin", args -> {if (args.size()==1) {return Math.sin(args.get(0));} return Double.NaN;});
            put("cos", args -> {if (args.size()==1) {return Math.cos(args.get(0));} return Double.NaN;});
            put("tan", args -> {if (args.size()==1) {return Math.tan(args.get(0));} return Double.NaN;});
            put("tg", args -> {if (args.size()==1) {return Math.tan(args.get(0));} return Double.NaN;});
            put("sec", args -> {if (args.size()==1) {return Functions.sec(args.get(0));} return Double.NaN;});
            put("csc", args -> {if (args.size()==1) {return Functions.csc(args.get(0));} return Double.NaN;});
            put("sinh", args -> {if (args.size()==1) {return Math.sinh(args.get(0));} return Double.NaN;});
            put("cosh", args -> {if (args.size()==1) {return Math.cosh(args.get(0));} return Double.NaN;});
            put("tanh", args -> {if (args.size()==1) {return Math.tanh(args.get(0));} return Double.NaN;});
            put("coth", args -> {if (args.size()==1) {return Functions.coth(args.get(0));} return Double.NaN;});
            put("sech", args -> {if (args.size()==1) {return Functions.sech(args.get(0));} return Double.NaN;});
            put("csch", args -> {if (args.size()==1) {return Functions.csch(args.get(0));} return Double.NaN;});
            put("sinp", args -> {if (args.size()==1) {return Functions.sinp(args.get(0));} return Double.NaN;});
            put("cosp", args -> {if (args.size()==1) {return Functions.cosp(args.get(0));} return Double.NaN;});
            put("tanp", args -> {if (args.size()==1) {return Functions.tanp(args.get(0));} return Double.NaN;});
            put("cotp", args -> {if (args.size()==1) {return Functions.cotp(args.get(0));} return Double.NaN;});
            put("secp", args -> {if (args.size()==1) {return Functions.secp(args.get(0));} return Double.NaN;});
            put("cscp", args -> {if (args.size()==1) {return Functions.cscp(args.get(0));} return Double.NaN;});
            put("arcsin", args -> {if (args.size()==1) {return Math.asin(args.get(0));} return Double.NaN;});
            put("arccos", args -> {if (args.size()==1) {return Math.acos(args.get(0));} return Double.NaN;});
            put("arctan", args -> {if (args.size()==1) {return Math.atan(args.get(0));} return Double.NaN;});
            put("arccot", args -> {if (args.size()==1) {return Functions.acot(args.get(0));} return Double.NaN;});
            put("arcsec", args -> {if (args.size()==1) {return Functions.asec(args.get(0));} return Double.NaN;});
            put("arccsc", args -> {if (args.size()==1) {return Functions.acsc(args.get(0));} return Double.NaN;});
            put("arcsinh", args -> {if (args.size()==1) {return Functions.asinh(args.get(0));} return Double.NaN;});
            put("arccosh", args -> {if (args.size()==1) {return Functions.acosh(args.get(0));} return Double.NaN;});
            put("arctanh", args -> {if (args.size()==1) {return Functions.atanh(args.get(0));} return Double.NaN;});
            put("arccoth", args -> {if (args.size()==1) {return Functions.acoth(args.get(0));} return Double.NaN;});
            put("arcsech", args -> {if (args.size()==1) {return Functions.asech(args.get(0));} return Double.NaN;});
            put("arccsch", args -> {if (args.size()==1) {return Functions.acsch(args.get(0));} return Double.NaN;});
            put("arcsinp", args -> {if (args.size()==1) {return Functions.asinp(args.get(0));} return Double.NaN;});
            put("arccosp", args -> {if (args.size()==1) {return Functions.acosp(args.get(0));} return Double.NaN;});
            put("arctanp", args -> {if (args.size()==1) {return Functions.atanp(args.get(0));} return Double.NaN;});
            put("arccotp", args -> {if (args.size()==1) {return Functions.acotp(args.get(0));} return Double.NaN;});
            put("arcsecp", args -> {if (args.size()==1) {return Functions.asecp(args.get(0));} return Double.NaN;});
            put("arccscp", args -> {if (args.size()==1) {return Functions.acscp(args.get(0));} return Double.NaN;});
            put("exp", args -> {if (args.size()==1) {return Math.exp(args.get(0));} return Double.NaN;});
            put("ln", args -> {if (args.size()==1) {return Math.log(args.get(0));} return Double.NaN;});
            put("log", args -> {if (args.size()==2) {return Functions.logb(args.get(0),args.get(1));} return Double.NaN;});
            put("pow", args -> {if (args.size()==2) {return Math.pow(args.get(0),args.get(1));} return Double.NaN;});
            put("Gamma", args -> {if (args.size()==1) {return Functions.gamma(args.get(0));} return Double.NaN;});
            put("fact", args -> {if (args.size()==1) {return Functions.factorial(args.get(0));} return Double.NaN;});
            put("factorial", args -> {if (args.size()==1) {return Functions.factorial(args.get(0));} return Double.NaN;});
            put("nCr", args -> {if (args.size()==2) {return Functions.nCr(args.get(0),args.get(1));} return Double.NaN;});
            put("nPr", args -> {if (args.size()==2) {return Functions.nPr(args.get(0),args.get(1));} return Double.NaN;});
            put("mean", Functions::mean);
            put("max", Functions::max);
            put("min", Functions::min);
            put("lcm", Functions::lcm);
            put("gcd", Functions::gcd);
            put("sign", args -> {if (args.size()==1) {return Functions.sgn(args.get(0));} return Double.NaN;});
            put("ceil", args -> {if (args.size()==1) {return Functions.ceil(args.get(0));} return Double.NaN;});
            put("floor", args -> {if (args.size()==1) {return Functions.floor(args.get(0));} return Double.NaN;});
            put("round", args -> {if (args.size()==1) {return Functions.round(args.get(0));} return Double.NaN;});
            put("psi", args -> {if (args.size()==1) {return Functions.diGamma(args.get(0));} return Double.NaN;});
            put("rect", args -> {if (args.size()==1) {return Functions.rect(args.get(0));} return Double.NaN;});
            put("Pi", args -> {if (args.size()==1) {return Functions.rect(args.get(0));} return Double.NaN;});
            put("tri", args -> {if (args.size()==1) {return Functions.tri(args.get(0));} return Double.NaN;});
            put("abs", args -> {if (args.size()==1) {return Math.abs(args.get(0));} return Double.NaN;});
            put("random", Functions::random);
            put("W", args -> {if (args.size()==1) {return Functions.W0(args.get(0));} return Double.NaN;});
            put("W_0", args -> {if (args.size()==1) {return Functions.W0(args.get(0));} return Double.NaN;});
            put("W_{0", args -> {if (args.size()==1) {return Functions.W0(args.get(0));} return Double.NaN;});
            put("SuSin", args -> {if (args.size()==1) {return Functions.SuSun(args.get(0));} return Double.NaN;});
            put("zeta", args -> {if (args.size()==1) {return Functions.zeta(args.get(0));} return Double.NaN;});
            put("tet_2", args -> {if (args.size()==1) {return Functions.tetBase2(args.get(0));} return Double.NaN;});
            put("tet_{2", args -> {if (args.size()==1) {return Functions.tetBase2(args.get(0));} return Double.NaN;});
            put("erf", args -> {if (args.size()==1) {return Functions.erf(args.get(0));} return Double.NaN;});
            put("invGamma", args -> {if (args.size()==1) {return Functions.invGamma(args.get(0));} return Double.NaN;});
            put("arcGamma", args -> {if (args.size()==1) {return Functions.invGamma(args.get(0));} return Double.NaN;});
            put("pi_{prime", args -> {if (args.size()==1) {return Functions.primePi(args.get(0));} return Double.NaN;});
            put("pi_{Prime", args -> {if (args.size()==1) {return Functions.primePi(args.get(0));} return Double.NaN;});
            put("mu", args -> {if (args.size()==1) {return Functions.muMobius(args.get(0));} return Double.NaN;});
            put("P", args -> {if (args.size()==1) {return Functions.primeZeta(args.get(0));} return Double.NaN;});
            put("sigma", args -> {
                if (args.size()==1) {return Functions.divisorSigma(args.get(0),1);}
                else if (args.size()==2) {return Functions.divisorSigma(args.get(0),args.get(1));}
                return Double.NaN;});
            put("isPrime", args -> {if (args.size()==1) {return Functions.isPrime(args.get(0));} return Double.NaN;});
            put("omega", args -> {if (args.size()==1) {return Functions.omegaPrime(args.get(0));} return Double.NaN;});
            put("Omega", args -> {if (args.size()==1) {return Functions.OmegaPrime(args.get(0));} return Double.NaN;});
            put("lambda", args -> {if (args.size()==1) {return Functions.lambdaLiouville(args.get(0));} return Double.NaN;});
            put("phi_{euler", args -> {if (args.size()==1) {return Functions.phiEuler(args.get(0));} return Double.NaN;});
            put("phi_{Euler", args -> {if (args.size()==1) {return Functions.phiEuler(args.get(0));} return Double.NaN;});
            put("collatzLen", args -> {if (args.size()==1) {return Functions.collatzLen(args.get(0));} return Double.NaN;});
            put("collatzMax", args -> {if (args.size()==1) {return Functions.collatzMax(args.get(0));} return Double.NaN;});
            put("C", args -> {if (args.size()==1) {return Functions.compositeZeta(args.get(0));} return Double.NaN;});
            put("xi", args -> {if (args.size()==1) {return Functions.xiRiemann(args.get(0));} return Double.NaN;});
            put("theta_{Chebyshev", args -> {if (args.size()==1) {return Functions.thetaChebyshev(args.get(0));} return Double.NaN;});
            put("Lambda", args -> {if (args.size()==1) {return Functions.lambdaMangoldt(args.get(0));} return Double.NaN;});
            put("psi_{Chebyshev", args -> {if (args.size()==1) {return Functions.psiChebyshev(args.get(0));} return Double.NaN;});
            put("B", args -> {
                if (args.size()==1) {return Functions.bernoulli(args.get(0));}
                else if (args.size()==2) {return Functions.beta(args.get(0),args.get(1));}
                return Double.NaN;});
            put("Beta", args -> {if (args.size()==2) {return Functions.beta(args.get(0),args.get(1));} return Double.NaN;});
            put("F", args -> {if (args.size()==1) {return Functions.fibonacci(args.get(0));} return Double.NaN;});
            put("Fibonacci", args -> {if (args.size()==1) {return Functions.fibonacci(args.get(0));} return Double.NaN;});
            put("ssqrt", args -> {if (args.size()==1) {return Functions.ssqrt(args.get(0));} return Double.NaN;});
            put("sqsrt", args -> {if (args.size()==1) {return Functions.ssqrt(args.get(0));} return Double.NaN;});
            put("trunc", args -> {if (args.size()==1) {return Functions.trunc(args.get(0));} return Double.NaN;});
            put("arcFact", args -> {if (args.size()==1) {return Functions.invFactorial(args.get(0));} return Double.NaN;});
            put("arcFactorial", args -> {if (args.size()==1) {return Functions.invFactorial(args.get(0));} return Double.NaN;});
            put("invFact", args -> {if (args.size()==1) {return Functions.invFactorial(args.get(0));} return Double.NaN;});
            put("invFactorial", args -> {if (args.size()==1) {return Functions.invFactorial(args.get(0));} return Double.NaN;});
            put("Fermi", args -> {
                if (args.size()==1) {return Functions.Fermi(args.get(0));}
                else if (args.size()==2) {return Functions.Fermi(args.get(0),args.get(1));}
                return Double.NaN;});
            put("dms", args -> {if (args.size()==1) {return Functions.dms(args.get(0));} return Double.NaN;});
            put("deg", args -> {if (args.size()==1) {return Functions.deg(args.get(0));} return Double.NaN;});
            //put("arg", args -> {if (args.size()==1) {return Functions.arg(args.get(0));} return Double.NaN;});
            put("Bernoulli", args -> {if (args.size()==1) {return Functions.bernoulli(args.get(0));} return Double.NaN;});
            put("bernoulli", args -> {if (args.size()==1) {return Functions.bernoulli(args.get(0));} return Double.NaN;});
            put("digits", args -> {
                if (args.size()==1) {return Functions.digits(args.get(0));}
                else if (args.size()==2) {return Functions.digits(args.get(0),args.get(1));}
                return Double.NaN;});
            put("concat", args -> {
                if (args.size()==2) {return Functions.concat(args.get(0),args.get(1));}
                else if (args.size()==3) {return Functions.concat(args.get(0),args.get(1),args.get(2));}
                return Double.NaN;});
            put("xexp", args -> {if (args.size()==1) {return Functions.xexp(args.get(0));} return Double.NaN;});
            put("expinf", args -> {if (args.size()==1) {return Functions.tetHeightInf(args.get(0));} return Double.NaN;});
            put("harmonic", args -> {if (args.size()==1) {return Functions.harmonic(args.get(0));} return Double.NaN;});
            put("Harmonic", args -> {if (args.size()==1) {return Functions.harmonic(args.get(0));} return Double.NaN;});
            put("H", args -> {if (args.size()==1) {return Functions.harmonic(args.get(0));} return Double.NaN;});
            put("clamp", args -> {if (args.size()==3) {return Functions.clamp(args.get(0),args.get(1),args.get(2));} return Double.NaN;});
            /*put("digitsum", args -> {
                if (args.size()==1) {return Functions.digitsum(args.get(0));}
                else if (args.size()==2) {return Functions.digitsum(args.get(0),args.get(1));}
                return Double.NaN;});*/
            put("AuSin", args -> {if (args.size()==1) {return Functions.AuSin(args.get(0));} return Double.NaN;});
            put("sin_{iter", args -> {if (args.size()==2) {return Functions.sinIterated(args.get(0),args.get(1));} return Double.NaN;});
            put("tet_E", args -> {if (args.size()==1) {return Functions.tetBaseE(args.get(0));} return Double.NaN;});
            put("tet_{E", args -> {if (args.size()==1) {return Functions.tetBaseE(args.get(0));} return Double.NaN;});
            put("tet_e", args -> {if (args.size()==1) {return Functions.tetBaseE(args.get(0));} return Double.NaN;});
            put("tet_{e", args -> {if (args.size()==1) {return Functions.tetBaseE(args.get(0));} return Double.NaN;});
            put("tet", args -> {if (args.size()==1) {return Functions.tetBaseE(args.get(0));} return Double.NaN;});
            put("pen_E", args -> {if (args.size()==1) {return Functions.penBaseE(args.get(0));} return Double.NaN;});
            put("pen_{E", args -> {if (args.size()==1) {return Functions.penBaseE(args.get(0));} return Double.NaN;});
            put("pen_e", args -> {if (args.size()==1) {return Functions.penBaseE(args.get(0));} return Double.NaN;});
            put("pen_{e", args -> {if (args.size()==1) {return Functions.penBaseE(args.get(0));} return Double.NaN;});
            put("pen", args -> {if (args.size()==1) {return Functions.penBaseE(args.get(0));} return Double.NaN;});
            put("slog_E", args -> {if (args.size()==1) {return Functions.slogBaseE(args.get(0));} return Double.NaN;});
            put("slog_{E", args -> {if (args.size()==1) {return Functions.slogBaseE(args.get(0));} return Double.NaN;});
            put("slog_e", args -> {if (args.size()==1) {return Functions.slogBaseE(args.get(0));} return Double.NaN;});
            put("slog_{e", args -> {if (args.size()==1) {return Functions.slogBaseE(args.get(0));} return Double.NaN;});
            put("slog", args -> {if (args.size()==1) {return Functions.slogBaseE(args.get(0));} return Double.NaN;});
            put("sln", args -> {if (args.size()==1) {return Functions.slogBaseE(args.get(0));} return Double.NaN;});
            put("exp_{iter", args -> {if (args.size()==2) {return Functions.expIterated(args.get(0),args.get(1));} return Double.NaN;});
            put("pen_2", args -> {if (args.size()==1) {return Functions.penBase2(args.get(0));} return Double.NaN;});
            put("pen_{2", args -> {if (args.size()==1) {return Functions.penBase2(args.get(0));} return Double.NaN;});
            //put("DeltaLn", args -> {if (args.size()==2) {return Functions.deltaLn(args.get(0),args.get(1));} return Double.NaN;});
            put("SuFact", args -> {if (args.size()==1) {return Functions.SuFact(args.get(0));} return Double.NaN;});
            put("AuFact", args -> {if (args.size()==1) {return Functions.AuFact(args.get(0));} return Double.NaN;});
            put("fact_{iter", args -> {if (args.size()==2) {return Functions.factIterated(args.get(0),args.get(1));} return Double.NaN;});
            put("apen_E", args -> {if (args.size()==1) {return Functions.arcPenBaseE(args.get(0));} return Double.NaN;});
            put("apen_{E", args -> {if (args.size()==1) {return Functions.arcPenBaseE(args.get(0));} return Double.NaN;});
            put("apen_e", args -> {if (args.size()==1) {return Functions.arcPenBaseE(args.get(0));} return Double.NaN;});
            put("apen_{e", args -> {if (args.size()==1) {return Functions.arcPenBaseE(args.get(0));} return Double.NaN;});
            put("apen", args -> {if (args.size()==1) {return Functions.arcPenBaseE(args.get(0));} return Double.NaN;});
            put("tet_{iter", args -> {if (args.size()==2) {return Functions.tetEIterated(args.get(0),args.get(1));} return Double.NaN;});
            put("sinc", args -> {if (args.size()==1) {return Functions.sinc(args.get(0));} return Double.NaN;});
            put("cosc", args -> {if (args.size()==1) {return Functions.cosc(args.get(0));} return Double.NaN;});
            put("tanc", args -> {if (args.size()==1) {return Functions.tanc(args.get(0));} return Double.NaN;});
            put("cotc", args -> {if (args.size()==1) {return Functions.cotc(args.get(0));} return Double.NaN;});
            put("secc", args -> {if (args.size()==1) {return Functions.secc(args.get(0));} return Double.NaN;});
            put("cscc", args -> {if (args.size()==1) {return Functions.cscc(args.get(0));} return Double.NaN;});
            put("sinhc", args -> {if (args.size()==1) {return Functions.sinhc(args.get(0));} return Double.NaN;});
            put("coshc", args -> {if (args.size()==1) {return Functions.coshc(args.get(0));} return Double.NaN;});
            put("tanhc", args -> {if (args.size()==1) {return Functions.tanhc(args.get(0));} return Double.NaN;});
            put("cothc", args -> {if (args.size()==1) {return Functions.cothc(args.get(0));} return Double.NaN;});
            put("sechc", args -> {if (args.size()==1) {return Functions.sechc(args.get(0));} return Double.NaN;});
            put("cschc", args -> {if (args.size()==1) {return Functions.cschc(args.get(0));} return Double.NaN;});
            put("sinpc", args -> {if (args.size()==1) {return Functions.sinpc(args.get(0));} return Double.NaN;});
            put("cospc", args -> {if (args.size()==1) {return Functions.cospc(args.get(0));} return Double.NaN;});
            put("tanpc", args -> {if (args.size()==1) {return Functions.tanpc(args.get(0));} return Double.NaN;});
            put("cotpc", args -> {if (args.size()==1) {return Functions.cotpc(args.get(0));} return Double.NaN;});
            put("secpc", args -> {if (args.size()==1) {return Functions.secpc(args.get(0));} return Double.NaN;});
            put("cscpc", args -> {if (args.size()==1) {return Functions.cscpc(args.get(0));} return Double.NaN;});
            put("erfc", args -> {if (args.size()==1) {return Functions.erfc(args.get(0));} return Double.NaN;});
            put("tet_{e^1/e", args -> {if (args.size()==1) {return Functions.tetee1(args.get(0));} return Double.NaN;});
            put("pen_{e^1/e", args -> {if (args.size()==1) {return Functions.penee1(args.get(0));} return Double.NaN;});
            put("SuXexp", args -> {if (args.size()==1) {return Functions.SuXexp(args.get(0));} return Double.NaN;});
            put("AuXexp", args -> {if (args.size()==1) {return Functions.AuXexp(args.get(0));} return Double.NaN;});
            put("xexp_{iter", args -> {if (args.size()==2) {return Functions.xexpEIterated(args.get(0),args.get(1));} return Double.NaN;});
            put("qrsqrt", args -> {if (args.size()==1) {return Functions.Q_rsqrt(args.get(0));} return Double.NaN;});
            put("sqrt", args -> {if (args.size()==1) {return Math.sqrt(args.get(0));} return Double.NaN;});
            put("cbrt", args -> {if (args.size()==1) {return Math.cbrt(args.get(0));} return Double.NaN;});
            put("Si", args -> {if (args.size()==1) {return Functions.Si(args.get(0));} return Double.NaN;});
            put("Ci", args -> {if (args.size()==1) {return Functions.Ci(args.get(0));} return Double.NaN;});
            //put("big", args -> {if (args.size()==1) {return MainDraw.big(args.get(0));} return Double.NaN;});
            //put("small", args -> {if (args.size()==1) {return MainDraw.small(args.get(0));} return Double.NaN;});
            put("tg", args -> {if (args.size()==1) {return Math.tan(args.get(0));} return Double.NaN;});
            put("ctg", args -> {if (args.size()==1) {return Functions.cot(args.get(0));} return Double.NaN;});
            put("cosec", args -> {if (args.size()==1) {return Functions.csc(args.get(0));} return Double.NaN;});
            put("arctg", args -> {if (args.size()==1) {return Math.atan(args.get(0));} return Double.NaN;});
            put("arcctg", args -> {if (args.size()==1) {return Functions.acot(args.get(0));} return Double.NaN;});
            put("arccosec", args -> {if (args.size()==1) {return Functions.acsc(args.get(0));} return Double.NaN;});
            put("sh", args -> {if (args.size()==1) {return Math.sinh(args.get(0));} return Double.NaN;});
            put("ch", args -> {if (args.size()==1) {return Math.cosh(args.get(0));} return Double.NaN;});
            put("th", args -> {if (args.size()==1) {return Math.tanh(args.get(0));} return Double.NaN;});
            put("cth", args -> {if (args.size()==1) {return Functions.coth(args.get(0));} return Double.NaN;});
            put("sch", args -> {if (args.size()==1) {return Functions.sech(args.get(0));} return Double.NaN;});
            put("arsh", args -> {if (args.size()==1) {return Functions.asinh(args.get(0));} return Double.NaN;});
            put("arch", args -> {if (args.size()==1) {return Functions.acosh(args.get(0));} return Double.NaN;});
            put("arth", args -> {if (args.size()==1) {return Functions.atanh(args.get(0));} return Double.NaN;});
            put("arcth", args -> {if (args.size()==1) {return Functions.acoth(args.get(0));} return Double.NaN;});
            put("arsch", args -> {if (args.size()==1) {return Functions.asech(args.get(0));} return Double.NaN;});
            put("arcsch", args -> {if (args.size()==1) {return Functions.acsch(args.get(0));} return Double.NaN;});
            this.putAll(constFunctionMap);
        }};

    }

    public static class Lexeme {
        LexemeType type;
        String value;

        public Lexeme(LexemeType type, String value) {
            this.type = type;
            this.value = value;
        }

        public Lexeme(LexemeType type, Character value) {
            this.type = type;
            this.value = value.toString();
        }

        @Override
        public String toString() {
            return "Lexeme{" +
                    "type=" + type +
                    ", value='" + value + '\'' +
                    '}';
        }
    }

    public static class LexemeBuffer {
        private int pos;

        public List<Lexeme> lexemes;

        public LexemeBuffer(List<Lexeme> lexemes) {
            this.lexemes = lexemes;
            this.pos = 0;
        }

        public Lexeme next() {
            return lexemes.get(pos++);
        }

        public void back() {
            pos--;
        }

        public int getPos() {
            return pos;
        }
    }

    public static List<Lexeme> lexAnalyze(String expText) {
        ArrayList<Lexeme> lexemes = new ArrayList<>();
        int pos = 0;
        while (pos < expText.length()) {
            char c = expText.charAt(pos);
            switch (c) {
                case '(':
                    lexemes.add(new Lexeme(LexemeType.LEFT_BRACKET, c));
                    pos++;
                    continue;
                case ')':
                    lexemes.add(new Lexeme(LexemeType.RIGHT_BRACKET, c));
                    pos++;
                    continue;
                case '+':
                    lexemes.add(new Lexeme(LexemeType.OP_PLUS, c));
                    pos++;
                    continue;
                case '-':
                    lexemes.add(new Lexeme(LexemeType.OP_MINUS, c));
                    pos++;
                    continue;
                case '*':
                    lexemes.add(new Lexeme(LexemeType.OP_MUL, c));
                    pos++;
                    continue;
                case '/':
                    lexemes.add(new Lexeme(LexemeType.OP_DIV, c));
                    pos++;
                    continue;
                case ',':
                    lexemes.add(new Lexeme(LexemeType.COMMA, c));
                    pos++;
                    continue;
                default:
                    if ((c <= '9' && c >= '0') || c == '.') {
                        StringBuilder sb = new StringBuilder();
                        boolean expPresent = false, isNumChar;
                        char prevC = c;
                        do {
                            sb.append(c);
                            pos++;
                            if (pos >= expText.length()) {
                                break;
                            }
                            prevC = c;
                            c = expText.charAt(pos);
                            isNumChar = (c <= '9' && c >= '0') || c == '.';
                            if ((c == 'E' || c == 'e') && !expPresent) {
                                expPresent = true;
                                isNumChar = true;
                            }
                            if (c == '-' && (prevC == 'E' || prevC == 'e')) isNumChar = true;
                        } while (isNumChar);
                        lexemes.add(new Lexeme(LexemeType.NUMBER, sb.toString()));
                    } else {
                        if (c != ' ') {
                            if (
                                    (c >= 'a' && c <= 'z') ||
                                            (c >= 'A' && c <= 'Z') ||
                                            (c == '\\')
                            ) {
                                StringBuilder sb = new StringBuilder();
                                do {
                                    if (c != ' ') {
                                        sb.append(c);
                                    }
                                    pos++;
                                    if (pos >= expText.length()) {
                                        break;
                                    }
                                    c = expText.charAt(pos);
                                } while (
                                        (c >= '0' && c <= '9') ||
                                                (c >= 'a' && c <= 'z') ||
                                                (c >= 'A' && c <= 'Z') ||
                                                (c == '\\') || (c == ' ') || (c == '_') ||
                                                (c == '^') ||
                                                (c == '{') || (c == '}')
                                );
                                String K = sb.toString();
                                if (K.charAt(K.length() - 1) == '}') {
                                    K = K.substring(0, K.length() - 1);
                                }
                                K = K.replace("\\operatorname{", "");
                                K = K.replace("\\mathrm{", "");
                                K = K.replace("\\", "");
                                if (functionMap.containsKey(K)) {
                                    lexemes.add(new Lexeme(LexemeType.NAME, K));
                                    if (constFunctions.contains(K)) {
                                        lexemes.add(new Lexeme(LexemeType.LEFT_BRACKET, "("));
                                        lexemes.add(new Lexeme(LexemeType.RIGHT_BRACKET, ")"));
                                    }
                                } else {
                                    throw new RuntimeException("Unexpected function: " + K);
                                }
                            }
                        } else {
                            pos++;
                        }
                    }
            }
        }
        lexemes.add(new Lexeme(LexemeType.EOF, ""));
        return lexemes;

    }

    public static double expr(LexemeBuffer lexemes) {
        Lexeme lexeme = lexemes.next();
        if (lexeme.type == LexemeType.EOF) {
            return 0;
        } else {
            lexemes.back();
            return plusminus(lexemes);
        }
    }

    public static double plusminus(LexemeBuffer lexemes) {
        double value = multdiv(lexemes);
        while (true) {
            Lexeme lexeme = lexemes.next();
            switch (lexeme.type) {
                case OP_PLUS:
                    value += multdiv(lexemes);
                    break;
                case OP_MINUS:
                    value -= multdiv(lexemes);
                    break;
                case EOF:
                case COMMA:
                case RIGHT_BRACKET:
                    lexemes.back();
                    return value;
                default:
                    throw new RuntimeException("Unexpected token: " + lexeme.value
                            + " at position: " + lexemes.getPos());
            }
        }
    }

    public static double multdiv(LexemeBuffer lexemes) {
        double value = factor(lexemes);
        while (true) {
            Lexeme lexeme = lexemes.next();
            switch (lexeme.type) {
                case OP_MUL:
                    value *= factor(lexemes);
                    break;
                case OP_DIV:
                    value /= factor(lexemes);
                    break;
                case EOF:
                case RIGHT_BRACKET:
                case OP_PLUS:
                case COMMA:
                case OP_MINUS:
                    lexemes.back();
                    return value;
                default:
                    throw new RuntimeException("Unexpected token: " + lexeme.value
                            + " at position: " + lexemes.getPos());
            }
        }
    }

    public static double factor(LexemeBuffer lexemes) {
        Lexeme lexeme = lexemes.next();
        switch (lexeme.type) {
            case NAME:
                lexemes.back();
                return func(lexemes);
            case OP_MINUS:
                return -factor(lexemes);
            case NUMBER:
                return Double.parseDouble(lexeme.value);
            case LEFT_BRACKET:
                double value = plusminus(lexemes);
                lexeme = lexemes.next();
                if (lexeme.type != LexemeType.RIGHT_BRACKET) {
                    throw new RuntimeException("Unexpected token: " + lexeme.value
                            + " at position: " + lexemes.getPos());
                }
                return value;
            default:
                throw new RuntimeException("Unexpected token: " + lexeme.value
                        + " at position: " + lexemes.getPos());
        }
    }

    public static double func(LexemeBuffer lexemeBuffer) {
        String name = lexemeBuffer.next().value;
        Lexeme lexeme = lexemeBuffer.next();

        if (lexeme.type != LexemeType.LEFT_BRACKET) {
            throw new RuntimeException("Wrong function call syntax at " + lexeme.value);
        }

        ArrayList<Double> args = new ArrayList<>();

        lexeme = lexemeBuffer.next();
        if (lexeme.type != LexemeType.RIGHT_BRACKET) {
            lexemeBuffer.back();
            do {
                args.add(expr(lexemeBuffer));
                lexeme = lexemeBuffer.next();

                if (lexeme.type != LexemeType.COMMA && lexeme.type != LexemeType.RIGHT_BRACKET) {
                    throw new RuntimeException("Wrong function call syntax at " + lexeme.value);
                }

            } while (lexeme.type == LexemeType.COMMA);
        }
        return functionMap.get(name).apply(args);
    }
}