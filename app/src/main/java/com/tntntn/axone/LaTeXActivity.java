package com.tntntn.axone;

import android.app.Activity;
import android.content.Intent;
import android.os.Bundle;
import android.text.Editable;
import android.text.TextWatcher;
import android.util.Base64;
import android.util.Log;
import android.webkit.WebSettings;
import android.webkit.WebView;
import android.webkit.WebViewClient;
import android.widget.EditText;

import androidx.core.graphics.Insets;
import androidx.core.view.ViewCompat;
import androidx.core.view.WindowInsetsCompat;
// import android.content.SharedPreferences;
// import android.content.SharedPreferences;
// import android.content.SharedPreferences;
// import android.content.SharedPreferences;
// import android.content.SharedPreferences;

public class LaTeXActivity extends Activity {

    WebView webView;
    private volatile boolean pageLoaded = false;
    private String latexFormula = "\\\\(\\\\)";

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        Intent intent = new Intent(this, MainDraw.class);

        setContentView(R.layout.activity_la_te_xactivity);
        webView = (WebView) findViewById(R.id.mathWebView);
        webView.setClickable(false);
        webView.setLongClickable(false);
        WebSettings webSettings = webView.getSettings();
        webSettings.setLoadWithOverviewMode(true);
        webSettings.setJavaScriptEnabled(true);
        webSettings.setUseWideViewPort(true);

        String encodedHtml = Base64.encodeToString(html.getBytes(), Base64.NO_PADDING);

        webView.loadData(encodedHtml, "text/html", "base64");
        //webView.loadUrl("file:///assets/MathTemplate.html");
        webView.setWebViewClient(new WebViewClient() {
            @Override
            public void onPageFinished(WebView view, String url) {
                pageLoaded = true;
                webView.loadUrl("javascript:showFormula('" + LaTeXActivity.this.latexFormula  + "')");
                super.onPageFinished(view, url);
            }
        });

        ViewCompat.setOnApplyWindowInsetsListener(findViewById(R.id.main), (v, insets) -> {
            Insets systemBars = insets.getInsets(WindowInsetsCompat.Type.systemBars());
            v.setPadding(systemBars.left, systemBars.top, systemBars.right, systemBars.bottom);
            return insets;
        });
        // TO DO: Control zoom, dx, dy with mouse/finger instead of input fields.
        final EditText latexText = findViewById(R.id.editLatexText);
        final EditText dxText = findViewById(R.id.dxInput);
        final EditText dyText = findViewById(R.id.dyInput);
        final EditText zoomText = findViewById(R.id.zoomInput);
        final EditText colorText = findViewById(R.id.colorInput);
        final EditText thickText = findViewById(R.id.thickInput);
        latexText.setText(MainDraw.formula);
        dxText.setText(MainDraw.dx1);
        dyText.setText(MainDraw.dy1);
        zoomText.setText(MainDraw.zoom1);
        colorText.setText(MainDraw.color1);
        thickText.setText(MainDraw.thick1);
        setFormula("`"+String.valueOf(latexText.getText()).replace("\\","\\\\")+"`");

        latexText.addTextChangedListener(new TextWatcher() {
            public void afterTextChanged(Editable s) {
                String latexTextText = String.valueOf(latexText.getText()).replace("\\","\\\\");
                setFormula("`"+latexTextText+"`");
                MainDraw.formula=latexTextText;
            }
            public void beforeTextChanged(CharSequence s, int start, int count, int after) {}
            public void onTextChanged(CharSequence s, int start, int before, int count) {} {}
        });

        dxText.addTextChangedListener(new TextWatcher() {
            public void afterTextChanged(Editable s) {
                MainDraw.dx1=String.valueOf(dxText.getText());
            }
            public void beforeTextChanged(CharSequence s, int start, int count, int after) {}
            public void onTextChanged(CharSequence s, int start, int before, int count) {} {}
        });

        dyText.addTextChangedListener(new TextWatcher() {
            public void afterTextChanged(Editable s) {
                MainDraw.dy1=String.valueOf(dyText.getText());
            }
            public void beforeTextChanged(CharSequence s, int start, int count, int after) {}
            public void onTextChanged(CharSequence s, int start, int before, int count) {} {}
        });

        zoomText.addTextChangedListener(new TextWatcher() {
            public void afterTextChanged(Editable s) {
                MainDraw.zoom1=String.valueOf(zoomText.getText());
            }
            public void beforeTextChanged(CharSequence s, int start, int count, int after) {}
            public void onTextChanged(CharSequence s, int start, int before, int count) {} {}
        });

        colorText.addTextChangedListener(new TextWatcher() {
            public void afterTextChanged(Editable s) {
                MainDraw.color1=String.valueOf(colorText.getText());
            }
            public void beforeTextChanged(CharSequence s, int start, int count, int after) {}
            public void onTextChanged(CharSequence s, int start, int before, int count) {} {}
        });

        thickText.addTextChangedListener(new TextWatcher() {
            public void afterTextChanged(Editable s) {
                MainDraw.thick1=String.valueOf(thickText.getText());
            }
            public void beforeTextChanged(CharSequence s, int start, int count, int after) {}
            public void onTextChanged(CharSequence s, int start, int before, int count) {} {}
        });
    }

    //public LaTeXActivity(int context) {
    //    super(context);
   // }

    public void setFormula(String latexFormula) {
        this.latexFormula = latexFormula;
        if (pageLoaded) {
            webView.loadUrl("javascript:showFormula('" + LaTeXActivity.this.latexFormula  + "')");
        } else {
            Log.e(LaTeXActivity.class.getSimpleName(), "Page is not loaded yet.");
        }
    }

    private static final String html = "<!DOCTYPE html>\n" +
            "<html>\n" +
            "<head>\n" +
            "    <meta charset=\"UTF-8\">\n" +
            "    <meta name=\"viewport\"\n" +
            "          content=\"\n" +
            "              width=device-width,\n" +
            "              height=device-height,\n" +
            "              initial-scale=1,\n" +
            "              user-scalable=0\" />\n" +
            "    <title>MathJax AsciiMath Template</title>\n" +
            "\n" +
            "    <script type=\"text/x-mathjax-config\">\n" +
            "  MathJax.Hub.Config({\n" +
            "  messageStyle: \"none\",\n" +
            "  /*styles: {\n" +
            "    \"#MathJax_Message\": {left: \"\", right: 0},\n" +
            "    \"#MathJax_MSIE_Frame\": {left: \"\", right: 0}\n" +
            "  },*/\n" +
            "  showProcessingMessages: false,\n" +
            "    tex2jax: { inlineMath: [['$','$'],['\\\\(','\\\\)']] },\n" +
            "    asciimath2jax: {\n" +
            "      delimiters: [['`','`'], ['$','$']]\n" +
            "    },\n" +
            "    \"fast-preview\": {\n" +
            "    Chunks: {EqnChunk: 10000, EqnChunkFactor: 1, EqnChunkDelay: 0},\n" +
            "    color: \"inherit!important\",\n" +
            "    updateTime: 30, updateDelay: 6,\n" +
            "    messageStyle: \"none\",\n" +
            "    disabled: false\n" +
            "  }\n" +
            "  });\n" +
            "</script>\n" +
            "\n" +
            "    <script type=\"text/javascript\"\n" +
            "            src=\"https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_HTMLorMML\"></script>\n" +
            "\n" +
            "    <script>\n" +
            "  (function () {\n" +
            "    var QUEUE = MathJax.Hub.queue;  // shorthand for the queue\n" +
            "    MathJax.Hub.processSectionDelay = 0;\n" +
            "  //  var math = null;                // the element jax for the math output.\n" +
            "\n" +
            "/*\n" +
            "    QUEUE.Push(function () {\n" +
            "      math = MathJax.Hub.getAllJax(\"ascii_output\")[0];\n" +
            "    });\n" +
            "*/\n" +
            "    window.UpdateMath = function (TeX) {\n" +
            "      QUEUE.Push([\"Typeset\", MathJax.Hub, \"ascii_output\"]);\n" +
            "    }\n" +
            "  })();\n" +
            "\n" +
            "  function showFormula (formula) {\n" +
            "    //MathJax.Message.File = console.log;\n" +
            "    document.getElementById(\"ascii_output\").textContent = formula;\n" +
            "    window.UpdateMath(formula);\n" +
            "  }\n" +
            "  </script>\n" +
            "</head>\n" +
            "<body style=\"text-align: center;\" onload=\"showFormula('\\\\(\\\\)')\">\n" +
            "<span id=\"ascii_output\">\n" +
            "</span>\n" +
            "</body>\n" +
            "</html>\n";

}