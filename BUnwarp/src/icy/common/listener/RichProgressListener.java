package icy.common.listener;

/**
 * Progress notification listener enriched with optional comments and optional data.
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public interface RichProgressListener {
	public boolean notifyProgress(double position, double length, String comment, Object data);
}
