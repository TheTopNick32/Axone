package com.tntntn.axone;

import android.content.Intent;
import android.os.Bundle;
import android.view.Menu;
import android.view.MenuItem;
import android.view.Window;

import androidx.appcompat.app.AppCompatActivity;

// TODO: Make a video for the project.
public class MainActivity extends AppCompatActivity {
    Object[] settings = {false, false, false, 0d, 0d, 1d};
    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        requestWindowFeature(Window.FEATURE_NO_TITLE);
        setContentView(new MainDraw(this,settings));
    }

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {

        getMenuInflater().inflate(R.menu.main_menu, menu);
        return true;
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        int id = item.getItemId();

        switch(id){
            case R.id.x_log_settings:
                settings[0]=!(boolean)(settings[0]);
                setContentView(new MainDraw(this, settings));
                return true;
            case R.id.y_log_settings:
                settings[1]=!(boolean)(settings[1]);
                setContentView(new MainDraw(this, settings));
                return true;
            // case R.id.polar_settings:
            //     settings[2]=!(boolean)(settings[2]);
            //     setContentView(new MainDraw(this, settings));
            //     return true;
            case R.id.input:
                Intent intent=new Intent(this,LaTeXActivity.class);
                startActivity(intent);
                return true;
            case R.id.tutorial:
                Intent intent2=new Intent(this, TutorialActivity.class);
                startActivity(intent2);
                return true;
        }
        return super.onOptionsItemSelected(item);
    }
}