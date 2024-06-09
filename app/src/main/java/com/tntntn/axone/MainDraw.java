package com.tntntn.axone;

import static java.lang.Integer.max;

import android.content.Context;
import android.graphics.Canvas;
import android.graphics.Color;
import android.graphics.Paint;
import android.view.View;

import java.util.List;

public class MainDraw extends View {
    private static String default_color = "#2fe91f";

    public static String formula = "";
    public static String dx1 = "0";
    public static String dy1 = "0";
    public static String zoom1 = "1";
    public static String color1 = default_color;
    public static String thick1 = "5";
    public static int type;
    public static boolean inversed;
    Object[] settings;

    /*[x_log, y_log, is_polar, dx, dy, zoom]*/
    public MainDraw(Context context, Object[] settings) {
        super(context);
        this.settings = new Object[settings.length];

        System.arraycopy(settings, 0, this.settings, 0, settings.length);
    }

    Paint paint = new Paint();
    boolean draw = true;
    List<LaTeXReader.Lexeme> latexPrepared;
    List<LaTeXReader.Lexeme> latexPreparedX;
    List<LaTeXReader.Lexeme> latexPreparedY;
    @Override
    protected void onDraw(Canvas canvas) {
        formula = formula.replace("\n", "");
        int equal = formula.indexOf('=');
        try {
            if (equal==-1) {
                latexPrepared = LaTeXReader.latexPrepare(formula);
                type = 1;
                inversed = false;
            } else if (formula.substring(0, equal).equals("y") && !formula.substring(equal).contains("y")) {
                latexPrepared = LaTeXReader.latexPrepare(formula.substring(equal+1));
                type = 1;
                inversed = false;
            } else if (formula.substring(equal + 1).equals("y") && !formula.substring(0,equal).contains("y")) {
                latexPrepared = LaTeXReader.latexPrepare(formula.substring(0,equal));
                type = 1;
                inversed = false;
            } else if (formula.substring(0, equal).equals("x") && !formula.substring(equal).contains("x")) {
                latexPrepared = LaTeXReader.latexPrepare(formula.substring(equal+1));
                type = 1;
                inversed = true;
            } else if (formula.substring(equal + 1).equals("x") && !formula.substring(0,equal).contains("x")) {
                latexPrepared = LaTeXReader.latexPrepare(formula.substring(0,equal));
                type = 1;
                inversed = true;
            } else {
                latexPreparedX = LaTeXReader.latexPrepare(formula.substring(0, equal));
                latexPreparedY = LaTeXReader.latexPrepare(formula.substring(equal + 1));
                type = 0;
                inversed = false;
            }
            double q = latexPrepared!=null?LaTeXReader.latexCalc(latexPrepared,0,0):Double.NaN;
            q = latexPreparedX!=null?LaTeXReader.latexCalc(latexPreparedX,0,0):Double.NaN;
            q = latexPreparedY!=null?LaTeXReader.latexCalc(latexPreparedY,0,0):Double.NaN;
        } catch (Exception e) {
            displayError(e, canvas);
            draw=false;
        }
        super.onDraw(canvas);
        if (draw) {
            int w = getWidth();
            int h = getHeight();

            /*Settings*/
            boolean x_log = (boolean) settings[0];
            boolean y_log = (boolean) settings[1];
            boolean is_polar = (boolean) settings[2];
            double dx = 0;
            double dy = 0;
            double zoom = 1;
            //String color = "#2fe91f";
            String color = color1;
            int thick = 5;
            try {
                dx = (double) Double.parseDouble(dx1);
            } catch (Exception ex) {
                System.err.println("Cannot convert \"" + dx1 + "\" to double!\n" + ex.toString());
            }
            try {
                dy = (double) Double.parseDouble(dy1);
            } catch (Exception ex) {
                System.err.println("Cannot convert \"" + dy1 + "\" to double!\n" + ex.toString());
            }
            try {
                zoom = (double) Double.parseDouble(zoom1);
            } catch (Exception ex) {
                System.err.println("Cannot convert \"" + zoom1 + "\" to double!\n" + ex.toString());
            }
            try {
                thick = (int) Math.round(Double.parseDouble(thick1));
            } catch (Exception ex) {
                System.err.println("Cannot convert \"" + thick1 + "\" to double!\n" + ex.toString());
            }
            if (color.isEmpty()) {
                System.err.println("MainDraw.color1 is an empty string!");
                color = default_color;
                color1 = default_color;
            }
            if (color.charAt(0)!='#') {
                color="#"+color;
            }
            if (color.length()!=7) {
                System.err.println("Cannot use \"" + color1 + "\" as a color!");
                color = default_color;
            } else {
                for (int C = 1; C<7;C++) {
                    char c = color1.charAt(C);
                    if (!((c>='A'&&c<='F') || (c>='a'&&c<='f') || Character.isDigit(c))) {
                        System.err.println("Cannot use \"" + color1 + "\" as a color!");
                        color = default_color;
                        break;
                    }
                }
            }
            thick = max(1,thick);
            zoom = zoom * w / 20;
            if (x_log || y_log) {
                is_polar = false;
            }

            /*Plane*/
            if (!is_polar) {
                drawPlane(w, h, dx, dy, zoom, x_log, y_log, canvas);
            } else {
                drawPolarPlane(w, h, dx, dy, zoom, canvas);
            }

            /*Graph f(y)=f(x)*/
            // TO DO: Very optimise
            if (type == 0 && !formula.isEmpty()) {
                paint.setColor(Color.parseColor(color));
                for (int i = -w / 2; i < w / 2; i++) {
                    for (int j = -h / 2; j < h / 2; j++) {
                        double[][] aN = {
                                f(dx, dy, zoom, i, j - 1, x_log, y_log),
                                f(dx, dy, zoom, i + 1, j - 1, x_log, y_log),
                                f(dx, dy, zoom, i, j, x_log, y_log),
                                f(dx, dy, zoom, i + 1, j, x_log, y_log),
                        };
                        boolean[] a = {aN[0][0] < aN[0][1], aN[1][0] < aN[1][1], aN[2][0] < aN[2][1], aN[3][0] < aN[3][1]};
                        if (!(Double.isNaN(aN[0][0]) && Double.isNaN(aN[1][0]) && Double.isNaN(aN[2][0]) && Double.isNaN(aN[3][0])) && !(a[0] && a[1] && a[2] && a[3]) && !(!a[0] && !a[1] && !a[2] && !a[3])) {
                            canvas.drawCircle(i + w / 2f, j + h / 2f, (float) thick, paint);
                        }
                    }
                }
            }

            /*Graph y=f(x)*/
            if (type == 1 && !formula.isEmpty()) {
                paint.setColor(Color.parseColor(color));
                paint.setAlpha(255);
                paint.setStrokeWidth((float) (2 * thick));
                double step = 1 / zoom;
                double k = inversed ? (double) w / (double) h : 1;
                for (double t = dx - w / (2 * zoom * k) + step; t <= dx + w / (2 * zoom * k); t += step) {
                    canvas.drawLine((float) ((f(t, x_log, y_log, inversed)[0] - dx) * zoom + w / 2d),
                            (float) ((dy - f(t, x_log, y_log, inversed)[1]) * zoom + h / 2d),
                            (float) ((f(t - step, x_log, y_log, inversed)[0] - dx) * zoom + w / 2d),
                            (float) ((dy - f(t - step, x_log, y_log, inversed)[1]) * zoom + h / 2d),
                            paint);
                }
            }

            /*Markers*/
            if (is_polar) {
                drawPolarMarkers(w, h, dx, dy, zoom, canvas);
            }
            drawMarkers(w, h, dx, dy, zoom, x_log, y_log, canvas);
            /*String latex = "\\begin{array}{|c|l|||r|c|}";
            latex += "\\hline";
            latex += "\\text{Matrix}&\\multicolumn{2}{|c|}{\\text{Multicolumns}}&\\text{Font sizes commands}\\cr";
            latex += "\\hline";
            latex += "\\begin{pmatrix}\\alpha_{11}&\\cdots&\\alpha_{1n}\\cr\\hdotsfor{3}\\cr\\alpha_{n1}&\\cdots&\\alpha_{nn}\\end{pmatrix}&\\Large \\text{Large Right}&\\small \\text{small Left}&\\tiny \\text{tiny Tiny}\\cr";
            latex += "\\hline";
            latex += "\\multicolumn{4}{|c|}{\\Huge \\text{Huge Multicolumns}}\\cr";
            latex += "\\hline";
            latex += "\\end{array}";
            androidx.compose.ui.graphics.Color a1=new androidx.compose.ui.graphics.Color();
            TeXFormula formula = new TeXFormula(latex);
            formula.createPNG(TeXConstants.STYLE_DISPLAY, 20, "target/Example5.png", a1, a1);*/
        }
    }

    public double[] f(double dx, double dy, double zoom, int x0, int y0, boolean xl, boolean yl) {
        double x = dx + (double) x0 / zoom;
        double y = dy + (double) -y0 / zoom;
        if (xl) {
            x = Math.pow(10, x);
        }
        if (yl) {
            y = Math.pow(10, y);
        }
        /*return new double[]{74*(x*x)+52*(y*y)+489*x+127*y-33*x*y, 65000};*/
        /*return new double[]{((x + 2) * (x + 2) + (y - 2) * (y - 2) - 3) * ((x - 2) * (x - 2) + (y - 2) * (y - 2) - 3) * (x * x + y * y - 25) * (x * x + y * y + 2.5 * y / Math.abs(y)) * ((x + 1.8) * (x + 1.8) + (y - 1.8) * (y - 1.8) - 2) * ((x - 1.8) * (x - 1.8) + (y - 1.8) * (y - 1.8) - 2), 0};*/
        /*return new double[]{Math.pow(x*x+y*y-1,3), Math.E*x*x*y*y*y};*/
        /*return new double[]{Math.tan(x*x+y*y), 0};*/
        return new double[]{LaTeXReader.latexCalc(latexPreparedX, x, y), LaTeXReader.latexCalc(latexPreparedY, x, y)};
    }

    public double[] f(double t, boolean xl, boolean yl, boolean inverse) {
        //return new double[]{t,LaTeXReader.latexCalc("sin(sin(3.14-0.57)) * arctan(log(10, 8)) - \n * \n".replace("\n",Double.toString(t))),-10,10};
        if (inverse) {
            double t1=yl?Math.pow(10, t):t;
            double k = LaTeXReader.latexCalc(latexPrepared, t, t1);
            k=xl?Math.log10(k):k;
            return new double[]{k, t};
        } else {
            double t1=xl?Math.pow(10, t):t;
            double k = LaTeXReader.latexCalc(latexPrepared, t1, t);
            k=yl?Math.log10(k):k;
            return new double[]{t, k};
        }
    }
    public void drawPlane(int w, int h, double dx, double dy, double zoom, boolean xl, boolean yl, Canvas canvas) {
        Paint paint = new Paint();
        int Nx = 10;
        int Ny = (10 * h) / w;
        double small = small(w / zoom);
        double big = big(w / zoom);
        float dsx = (float) dx % (float) big;
        float dsy = (float) dy % (float) big;

        /*Small grid*/
        paint.setStrokeWidth(2);
        paint.setColor(Color.parseColor("#c0c0c0"));
        for (int i = 0; i < 5 * Nx; i++) {
            canvas.drawLine(w / 2f - (float) (i * small + dsx) * (float) zoom, 0, w / 2f - (float) (i * small + dsx) * (float) zoom, h, paint);
        }
        for (int i = 0; i < 5 * Nx; i++) {
            canvas.drawLine(w / 2f + (float) (i * small - dsx) * (float) zoom, 0, w / 2f + (float) (i * small - dsx) * (float) zoom, h, paint);
        }
        for (int i = 0; i < 5 * Ny; i++) {
            canvas.drawLine(0, h / 2f + (float) (i * small + dsy) * (float) zoom, w, h / 2f + (float) (i * small + dsy) * (float) zoom, paint);
        }
        for (int i = 0; i < 5 * Ny; i++) {
            canvas.drawLine(0, h / 2f - (float) (i * small - dsy) * (float) zoom, w, h / 2f - (float) (i * small - dsy) * (float) zoom, paint);
        }

        /*Big grid*/
        paint.setStrokeWidth(4);
        paint.setColor(Color.parseColor("#606060"));
        for (int i = 0; i < Nx; i++) {
            canvas.drawLine(w / 2f - (float) (i * big + dsx) * (float) zoom, 0, w / 2f - (float) (i * big + dsx) * (float) zoom, h, paint);
        }
        for (int i = 0; i < Nx; i++) {
            canvas.drawLine(w / 2f + (float) (i * big - dsx) * (float) zoom, 0, w / 2f + (float) (i * big - dsx) * (float) zoom, h, paint);
        }
        for (int i = 0; i < Ny; i++) {
            canvas.drawLine(0, h / 2f + (float) (i * big + dsy) * (float) zoom, w, h / 2f + (float) (i * big + dsy) * (float) zoom, paint);
        }
        for (int i = 0; i < Ny; i++) {
            canvas.drawLine(0, h / 2f - (float) (i * big - dsy) * (float) zoom, w, h / 2f - (float) (i * big - dsy) * (float) zoom, paint);
        }

        /*Axes*/
        paint.setStrokeWidth(8);
        paint.setColor(Color.parseColor("#000000"));
        if (!xl) {
            canvas.drawLine(w / 2f - (float) dx * (float) zoom, 0, w / 2f - (float) dx * (float) zoom, h, paint);
        }
        if (!yl) {
            canvas.drawLine(0, h / 2f + (float) dy * (float) zoom, w, h / 2f + (float) dy * (float) zoom, paint);
        }
    }

    // TODO: Show marks even if dy is "big".
    public void drawMarkers(int w, int h, double dx, double dy, double zoom, boolean xl, boolean yl, Canvas canvas) {
        Paint paint = new Paint();
        paint.setTextSize(w / 30f);
        int Nx = 10;
        int Ny = (10 * h) / w;
        double big = big(w / zoom);
        int k = (int) Math.floor(Math.abs(dx/big));
        int l = (int) Math.floor(Math.abs(dy/big));
        float xmyp = yl?h-20:h / 2f + 45 + (float) dy * (float) zoom;
        boolean y_out = (w / 2f - 10 - (float) dx * (float) zoom<15)||(w / 2f - 10 - (float) dx * (float) zoom>w-20);
        xmyp=(float) Functions.clamp(xmyp,30,h-15);
        float y_out_val = (y_out&&dx<0)?w:0;

        /*Text outline*/
        paint.setStyle(Paint.Style.STROKE);
        paint.setStrokeWidth(8);
        paint.setColor(Color.WHITE);
        paint.setTextAlign(Paint.Align.CENTER);
        for (int i =  max(0,k-10); i < k+Nx; i++) {
            if (i * big == 0) {
                continue;
            }
            canvas.drawText(doubleToString(xl?Math.pow(10,-i * big):-i * big), w / 2f - (float) (i * big + dx) * (float) zoom, xmyp, paint);
        }
        for (int i = max(0,k-10); i < k+Nx; i++) {
            if (i * big == 0) {
                continue;
            }
            canvas.drawText(doubleToString(xl?Math.pow(10,i * big):i * big), w / 2f - (float) (-i * big + dx) * (float) zoom, xmyp, paint);
        }
        if (!xl && !y_out) {paint.setTextAlign(Paint.Align.RIGHT);}
        else if (!xl && y_out_val>0) {paint.setTextAlign(Paint.Align.RIGHT);}
        else {paint.setTextAlign(Paint.Align.LEFT);}
        for (int i = max(0,l-10); i < l+Ny; i++) {
            if (i * big == 0) {
                continue;
            }
            canvas.drawText(doubleToString(yl?Math.pow(10,-i * big):-i * big), (xl||y_out)?y_out_val:w / 2f - 10 - (float) dx * (float) zoom, h / 2f + (float) (i * big) * (float) zoom + 15 + (float) dy * (float) zoom, paint);
        }
        for (int i = max(0,l-10); i < l+Ny; i++) {
            if (i * big == 0) {
                continue;
            }
            canvas.drawText(doubleToString(yl?Math.pow(10,i * big):i * big), (xl||y_out)?y_out_val:w / 2f - 10 - (float) dx * (float) zoom, h / 2f + (float) (-i * big) * (float) zoom + 15 + (float) dy * (float) zoom, paint);
        }
        if (!(xl&&yl)) {
            canvas.drawText((xl || yl) ? "1" : "0", w / 2f - 10 - (float) dx * (float) zoom, h / 2f + (float) dy * (float) zoom + 45, paint);
        }
        if (xl&&yl) {
            canvas.drawText("1", w / 2f - 10 - (float) dx * (float) zoom, xmyp, paint);
            canvas.drawText("1", 0, h / 2f + (float) dy * (float) zoom + 15, paint);
        }

        /*Text*/
        paint.setStyle(Paint.Style.FILL);
        paint.setStrokeWidth(0);
        paint.setColor(Color.BLACK);
        paint.setTextAlign(Paint.Align.CENTER);
        for (int i = max(0,k-10); i < k+Nx; i++) {
            if (i * big == 0) {
                continue;
            }
            canvas.drawText(doubleToString(xl?Math.pow(10,-i * big):-i * big), w / 2f - (float) (i * big + dx) * (float) zoom, xmyp, paint);
        }
        for (int i = max(0,k-10); i < k+Nx; i++) {
            if (i * big == 0) {
                continue;
            }
            canvas.drawText(doubleToString(xl?Math.pow(10,i * big):i * big), w / 2f - (float) (-i * big + dx) * (float) zoom, xmyp, paint);
        }
        if (!xl && !y_out) {paint.setTextAlign(Paint.Align.RIGHT);}
        else if (!xl && y_out_val>0) {paint.setTextAlign(Paint.Align.RIGHT);}
        else {paint.setTextAlign(Paint.Align.LEFT);}
        for (int i = max(0,l-10); i < l+Ny; i++) {
            if (i * big == 0) {
                continue;
            }
            canvas.drawText(doubleToString(yl?Math.pow(10,-i * big):-i * big), (xl||y_out)?y_out_val:w / 2f - 10 - (float) dx * (float) zoom, h / 2f + (float) (i * big) * (float) zoom + 15 + (float) dy * (float) zoom, paint);
        }
        for (int i = max(0,l-10); i < l+Ny; i++) {
            if (i * big == 0) {
                continue;
            }
            canvas.drawText(doubleToString(yl?Math.pow(10,i * big):i * big), (xl||y_out)?y_out_val:w / 2f - 10 - (float) dx * (float) zoom, h / 2f + (float) (-i * big) * (float) zoom + 15 + (float) dy * (float) zoom, paint);
        }
        if (!(xl&&yl)) {
            canvas.drawText((xl || yl) ? "1" : "0", w / 2f - 10 - (float) dx * (float) zoom, h / 2f + (float) dy * (float) zoom + 45, paint);
        }
        if (xl&&yl) {
            canvas.drawText("1", w / 2f - 10 - (float) dx * (float) zoom, xmyp, paint);
            canvas.drawText("1", 0, h / 2f + (float) dy * (float) zoom + 15, paint);
        }
    }

    public void drawPolarPlane(int w, int h, double dx, double dy, double zoom, Canvas canvas) {
        Paint paint = new Paint();
        double big = big(w / zoom);
        int Ny = (20 * h) / w;
        paint.setStyle(Paint.Style.STROKE);

        /*Radius*/
        paint.setStrokeWidth(4);
        paint.setColor(Color.parseColor("#606060"));
        for (int i = 0; i < Ny; i++) {
            canvas.drawCircle((float) (w / 2 - dx * zoom), (float) (h / 2 + dy * zoom), i * (float) big * (float) zoom, paint);
        }

        /*Angle*/
        paint.setStrokeWidth(2);
        paint.setColor(Color.parseColor("#c0c0c0"));
        for (double i = 0; i < Math.toRadians(360); i += Math.toRadians(15)) {
            canvas.drawLine((float) (w / 2 - dx * zoom), (float) (h / 2 + dy * zoom), (float) Math.cos(i) * (float) Math.max(2 * w, 2 * h) + (float) (w / 2 - dx * zoom), (float) Math.sin(i) * (float) Math.max(2 * w, 2 * h) + (float) (h / 2 + dy * zoom), paint);
        }

        /*Axes*/
        paint.setStrokeWidth(8);
        paint.setColor(Color.parseColor("#000000"));
        canvas.drawLine(w / 2f - (float) dx * (float) zoom, 0, w / 2f - (float) dx * (float) zoom, h, paint);
        canvas.drawLine(0, h / 2f + (float) dy * (float) zoom, w, h / 2f + (float) dy * (float) zoom, paint);
    }

    public void drawPolarMarkers(int w, int h, double dx, double dy, double zoom, Canvas canvas) {
        Paint paint = new Paint();
        paint.setTextSize(w / 25f);
        paint.setTextAlign(Paint.Align.CENTER);
        int N = 12;
        double big = big(w / zoom);
        float dsx = (float) dx * (float) zoom;
        float dsy = (float) dy * (float) zoom;
        float r = (float) (big * zoom * Functions.floor(w / (big * zoom * 2))) - 45;

        /*Text outline*/
        paint.setStyle(Paint.Style.STROKE);
        paint.setStrokeWidth(8);
        paint.setColor(Color.WHITE);
        for (double i = 0; i < Math.toRadians(360); i += Math.toRadians(360f / N)) {
            canvas.drawText(i == 0 ? "0" : one(ratioPi1(i), 1) + "π" + one(ratioPi2(i), 2), r * (float) Math.cos(i) + w / 2f - dsx, -r * (float) Math.sin(i) + h / 2f + 15 + dsy, paint);
        }

        /*Text*/
        paint.setStyle(Paint.Style.FILL);
        paint.setStrokeWidth(0);
        paint.setColor(Color.BLACK);
        for (double i = 0; i < Math.toRadians(360); i += Math.toRadians(360f / N)) {
            canvas.drawText(i == 0 ? "0" : one(ratioPi1(i), 1) + "π" + one(ratioPi2(i), 2), r * (float) Math.cos(i) + w / 2f - dsx, -r * (float) Math.sin(i) + h / 2f + 15 + dsy, paint);
        }

    }

    public static double big(double x) {
        double ans = 0;
        double x1 = Math.pow(10,Functions.frac(Math.log10(x)));
        if (1 <= x1 && x1 < 2) {
            ans = 0.2;
        }
        if (2 <= x1 && x1 < 5) {
            ans = 0.5;
        }
        if (5 <= x1 && x1 < 10) {
            ans = 1;
        }
        return ans * Math.pow(10, Functions.floor(Math.log10(x)));
    }

    public static double small(double x) {
        double ans = 0;
        double x1 = Math.pow(10,Functions.frac(Math.log10(x)));
        if (1 <= x1 && x1 < 2) {
            ans = 0.05;
        }
        if (2 <= x1 && x1 < 5) {
            ans = 0.1;
        }
        if (5 <= x1 && x1 < 10) {
            ans = 0.2;
        }
        return ans * Math.pow(10, Functions.floor(Math.log10(x)));
    }

    public static int ratioPi1(double x) {
        return (int) (Math.round(360 * x / Math.PI) / Functions.gcd(360, Math.round(360 * x / Math.PI)));
    }

    public static int ratioPi2(double x) {
        return (int) (360 / Functions.gcd(360, Math.round(360 * x / Math.PI)));
    }

    public static String one(int n, int p) {
        return p == 1 ? (n == 1 ? "" : Integer.toString(n)) : (n == 1 ? "" : "/" + n);
    }

    public static int[] pow_convert(String s, int index, boolean check_division) {
        int pow_start = 0;
        int pow_end = s.length()-1;
        char currchar;
        if (s.charAt(index-1)==')') {
            int counter=1;
            for (int i=index-2;i>=0;i--) {
                currchar=s.charAt(i);
                if (counter>0) {
                    if (currchar==')') {counter++;}
                    else if (currchar=='(') {counter--;}
                } else if (currchar==',' || currchar=='+' || currchar=='*' || currchar=='-' || currchar=='(') {
                    pow_start=i+1;
                    break;
                } else if (currchar=='/' && check_division) {
                    pow_start=pow_convert(s,index,false)[0];
                }
            }
        } else {
            for (int i=index;i>=0;i--) {
                currchar=s.charAt(i);
                if (currchar==',' || currchar=='+' || currchar=='*' || currchar=='-' || currchar=='(') {
                    pow_start=i+1;
                    break;
                } else if (currchar=='/') {
                    pow_start=pow_convert(s,index,false)[0];
                }
            }
        }
        if (Character.isDigit(s.charAt(index+1))) {
            pow_end=index+1;
        } else if (s.charAt(index+1)=='{') {
            int counter=1;
            for (int i=index+2;i<s.length();i++) {
                currchar=s.charAt(i);
                if (currchar=='{') {counter++;}
                else if (currchar=='}') {counter--;}
                pow_end=i-1;
                if (counter==0) {
                    break;
                }
            }
        } else {
            int counter=1;
            for (int i=index+2;i<s.length();i++) {
                currchar=s.charAt(i);
                if (currchar==',' || currchar=='+' || currchar=='*' || currchar=='-' || currchar=='/' || currchar=='(' || currchar==')') {
                    pow_end=i-1;
                    break;
                }
            }
        }
        return new int[] {pow_start,pow_end};
    }

    public static String pow_converter(String s, int start) {
        for (int i = start; i < s.length(); i++) {
            if (s.charAt(i)=='^') {
                int[] power = pow_convert(s,i,true);
                if (s.charAt(i+1)!='{') {return pow_converter(s.substring(0, power[0]) + "pow(" + s.substring(power[0], i) + "," + s.substring(i + 1, power[1] + 1) + ")" + s.substring(power[1] + 1), i);}
                else {return pow_converter(s.substring(0, power[0]) + "pow(" + s.substring(power[0], i) + "," + s.substring(i + 2, power[1] + 1) + ")" + s.substring(power[1] + 2), i);}
            }
        }
        return s;
    }

    public static String doubleToString(double a) {
        String k = Double.toString(a);
        if (Math.abs(a)==0) {
            return "0";
        }
        int v = k.indexOf('.');
        int v2 = k.indexOf('E');
        int i = v;
        if (v2!=-1) {
            while (i <= k.length() - 1) {
                if (k.charAt(i) == '0') {
                    k = k.substring(0, i) + k.substring(v2);
                    break;
                } else {
                    i++;
                }
            }
        }
        else {
            i++;
            boolean m=true;
            while (i <= k.length() - 1) {
                if (k.charAt(i) != '0' && m) {
                    m=false;
                    i++;
                } else if (k.charAt(i) == '0' && !m) {
                    k=(i<k.length()-2?k.substring(0, i+1):k);
                    break;
                } else {
                    i++;
                }
            }
        }
        if (k.charAt(k.length()-1)=='.') {
            k=k.substring(0,k.length()-1);
        }
        String STR = k.replace(".E","E");
        v2 = k.indexOf('E');
        if (v2==-1 && STR.length()>11) {
            STR=STR.substring(0,10);
        } else if (v2>11) {
            STR=STR.substring(0,10)+STR.substring(v2);
        }
        return STR;
    }

    public void displayError(Exception e, Canvas canvas) {
        Paint paint = new Paint();
        float size = getWidth() / 15f;
        paint.setTextSize(size);
        paint.setStyle(Paint.Style.FILL);
        paint.setStrokeWidth(4);
        paint.setColor(Color.RED);
        paint.setTextAlign(Paint.Align.LEFT);
        String error = e.toString();
        int to = (error.length()-1)/16;
        for (int i = 0; i < to; i++) {
            canvas.drawText(error.substring(i*16,(i+1)*16),size/4f,size*(i+1)*1.4f,paint);
        }
        canvas.drawText(error.substring(to*16),size/4f,size*(to+1)*1.4f,paint);
        System.err.println(e);
    }
}
