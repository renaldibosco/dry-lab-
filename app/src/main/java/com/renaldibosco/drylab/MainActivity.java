package com.renaldibosco.drylab;

import android.app.Activity;
import android.content.ActivityNotFoundException;
import android.content.ClipData;
import android.content.ClipboardManager;
import android.content.ContentValues;
import android.content.Context;
import android.content.Intent;
import android.graphics.Insets;
import android.net.Uri;
import android.os.Build;
import android.os.Bundle;
import android.provider.MediaStore;
import android.view.WindowInsets;
import android.webkit.JavascriptInterface;
import android.webkit.ValueCallback;
import android.webkit.WebChromeClient;
import android.webkit.WebResourceRequest;
import android.webkit.WebSettings;
import android.webkit.WebView;
import android.webkit.WebViewClient;
import android.widget.FrameLayout;
import android.widget.Toast;

/**
 * Dry Lab: the whole app (screens and biology engine) lives in assets/www and runs offline.
 * This class only shows it full screen and connects phone features: opening files, sharing, copying.
 */
public class MainActivity extends Activity {
    private static final int PICK_FILE = 42;
    private static final String HOME = "file:///android_asset/www/index.html";

    private WebView web;
    private ValueCallback<Uri[]> fileCallback;
    private Uri cameraUri;   // where the camera saves a new gel photo

    @Override
    protected void onCreate(Bundle state) {
        super.onCreate(state);

        FrameLayout root = new FrameLayout(this);
        root.setBackgroundColor(getColor(R.color.bar));
        web = new WebView(this);
        web.setBackgroundColor(getColor(R.color.page));
        root.addView(web, new FrameLayout.LayoutParams(
                FrameLayout.LayoutParams.MATCH_PARENT, FrameLayout.LayoutParams.MATCH_PARENT));
        setContentView(root);

        // Keep the page clear of the status bar, navigation bar and keyboard (Android 15 draws edge to edge).
        root.setOnApplyWindowInsetsListener((v, insets) -> {
            if (Build.VERSION.SDK_INT >= 30) {
                Insets b = insets.getInsets(WindowInsets.Type.systemBars() | WindowInsets.Type.ime());
                v.setPadding(b.left, b.top, b.right, b.bottom);
            } else {
                v.setPadding(insets.getSystemWindowInsetLeft(), insets.getSystemWindowInsetTop(),
                        insets.getSystemWindowInsetRight(), insets.getSystemWindowInsetBottom());
            }
            return insets;
        });

        WebSettings s = web.getSettings();
        s.setJavaScriptEnabled(true);
        s.setDomStorageEnabled(true);      // remembers your last sequence
        s.setAllowContentAccess(true);     // read files picked with the file chooser
        s.setAllowFileAccess(false);
        s.setDatabaseEnabled(true);
        web.addJavascriptInterface(new Bridge(), "Android");

        web.setWebViewClient(new WebViewClient() {
            @Override
            public boolean shouldOverrideUrlLoading(WebView view, WebResourceRequest request) {
                Uri u = request.getUrl();
                if ("file".equals(u.getScheme())) return false;
                try { startActivity(new Intent(Intent.ACTION_VIEW, u)); } catch (ActivityNotFoundException ignored) { }
                return true;
            }
        });

        web.setWebChromeClient(new WebChromeClient() {
            @Override
            public boolean onShowFileChooser(WebView view, ValueCallback<Uri[]> callback, FileChooserParams params) {
                if (fileCallback != null) fileCallback.onReceiveValue(null);
                fileCallback = callback;
                boolean image = false;
                String[] types = params.getAcceptTypes();
                if (types != null) for (String t : types) if (t != null && t.startsWith("image")) image = true;

                Intent pick = new Intent(Intent.ACTION_GET_CONTENT);
                pick.addCategory(Intent.CATEGORY_OPENABLE);
                pick.setType(image ? "image/*" : "*/*");   // .fasta and .gb files have no standard type
                Intent chooser = Intent.createChooser(pick, getString(image ? R.string.open_photo : R.string.open_file));

                // For gel photos also offer the camera (Android 10+: the photo goes to Pictures/Dry Lab, no permission needed)
                cameraUri = null;
                if (image && Build.VERSION.SDK_INT >= 29) {
                    try {
                        ContentValues v = new ContentValues();
                        v.put(MediaStore.Images.Media.DISPLAY_NAME, "drylab_gel_" + System.currentTimeMillis() + ".jpg");
                        v.put(MediaStore.Images.Media.MIME_TYPE, "image/jpeg");
                        v.put(MediaStore.Images.Media.RELATIVE_PATH, "Pictures/Dry Lab");
                        cameraUri = getContentResolver().insert(MediaStore.Images.Media.EXTERNAL_CONTENT_URI, v);
                        if (cameraUri != null) {
                            Intent cam = new Intent(MediaStore.ACTION_IMAGE_CAPTURE);
                            cam.putExtra(MediaStore.EXTRA_OUTPUT, cameraUri);
                            cam.addFlags(Intent.FLAG_GRANT_WRITE_URI_PERMISSION | Intent.FLAG_GRANT_READ_URI_PERMISSION);
                            if (cam.resolveActivity(getPackageManager()) != null) chooser.putExtra(Intent.EXTRA_INITIAL_INTENTS, new Intent[]{cam});
                            else { getContentResolver().delete(cameraUri, null, null); cameraUri = null; }
                        }
                    } catch (Exception e) { cameraUri = null; }
                }
                try {
                    startActivityForResult(chooser, PICK_FILE);
                } catch (ActivityNotFoundException e) {
                    fileCallback = null;
                    Toast.makeText(MainActivity.this, "No file manager found", Toast.LENGTH_SHORT).show();
                    return false;
                }
                return true;
            }
        });

        if (state != null) web.restoreState(state);
        if (web.getUrl() == null) web.loadUrl(HOME);
    }

    @Override
    protected void onActivityResult(int requestCode, int resultCode, Intent data) {
        if (requestCode == PICK_FILE && fileCallback != null) {
            Uri[] result = null;
            if (resultCode == RESULT_OK) {
                if (data != null && data.getData() != null) result = new Uri[]{data.getData()};
                else if (cameraUri != null) result = new Uri[]{cameraUri};          // photo just taken
            }
            if (result == null || (cameraUri != null && !cameraUri.equals(result[0]))) {
                if (cameraUri != null) try { getContentResolver().delete(cameraUri, null, null); } catch (Exception ignored) { }
            }
            cameraUri = null;
            fileCallback.onReceiveValue(result);
            fileCallback = null;
            return;
        }
        super.onActivityResult(requestCode, resultCode, data);
    }

    @Override
    protected void onSaveInstanceState(Bundle out) {
        super.onSaveInstanceState(out);
        web.saveState(out);
    }

    /** Back button: go back to the first tab before leaving the app. */
    @Override
    @SuppressWarnings("deprecation")
    public void onBackPressed() {
        web.evaluateJavascript("window.dryLabBack ? window.dryLabBack() : false", value -> {
            if (!"true".equals(value)) finish();
        });
    }

    @Override
    protected void onDestroy() {
        if (web != null) web.destroy();
        super.onDestroy();
    }

    /** Called from the page as window.Android.share(...) and window.Android.copy(...). */
    private class Bridge {
        @JavascriptInterface
        public void share(String title, String text) {
            Intent send = new Intent(Intent.ACTION_SEND);
            send.setType("text/plain");
            send.putExtra(Intent.EXTRA_SUBJECT, title);
            send.putExtra(Intent.EXTRA_TEXT, text);
            startActivity(Intent.createChooser(send, title));
        }

        @JavascriptInterface
        public void copy(String text) {
            runOnUiThread(() -> {
                ClipboardManager cm = (ClipboardManager) getSystemService(Context.CLIPBOARD_SERVICE);
                if (cm != null) cm.setPrimaryClip(ClipData.newPlainText("Dry Lab", text));
            });
        }
    }
}
