package main.java.smoge.utils;

import java.util.concurrent.TimeUnit;

public class Timer {
    private long start;

    public Timer() {
        this.start = System.currentTimeMillis();
    }

    public String getElapsedTime() {
        return timeFormatter(System.currentTimeMillis() - this.start);
    }

    private String timeFormatter(long millis) {
        return String.format("%02d:%02d:%02d", TimeUnit.MILLISECONDS.toHours(millis),
                TimeUnit.MILLISECONDS.toMinutes(millis) % TimeUnit.HOURS.toMinutes(1),
                TimeUnit.MILLISECONDS.toSeconds(millis) % TimeUnit.MINUTES.toSeconds(1));
    }
}
