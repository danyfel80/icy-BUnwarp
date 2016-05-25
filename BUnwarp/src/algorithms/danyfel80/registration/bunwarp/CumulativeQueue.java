package algorithms.danyfel80.registration.bunwarp;

import java.util.Vector;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class CumulativeQueue extends Vector<Double> {

	/**
	 * Generated Serial Version UID
	 */
	private static final long serialVersionUID = -3975876484175592536L;

	/** front index of the queue */
	private int ridx;
	/** rear index of the queue */
	private int widx;
	/** current length of the queue */
	private int currentLength;
	/** queue sum */
	private double sum;

	/*------------------------------------------------------------------*/
	/**
	 * Create a new instance of CumulativeQueue.
	 *
	 * @param length
	 *          length of the queue to be created
	 */
	public CumulativeQueue(int length) {
		currentLength = ridx = widx = 0;
		setSize(length);
	}

	/*------------------------------------------------------------------*/
	/**
	 * Get the current size of the queue.
	 *
	 * @return current size
	 */
	public int currentSize() {
		return currentLength;
	}

	/*------------------------------------------------------------------*/
	/**
	 * Get the sum of the queue.
	 *
	 * @return sum
	 */
	public double getSum() {
		return sum;
	}

	/*------------------------------------------------------------------*/
	/**
	 * Pop the value from the front of the queue.
	 *
	 * @return front value
	 */
	public double pop_front() {
		if (currentLength == 0)
			return 0.0;
		double x = ((Double) elementAt(ridx)).doubleValue();
		currentLength--;
		sum -= x;
		ridx++;
		if (ridx == size())
			ridx = 0;
		return x;
	}

	/*------------------------------------------------------------------*/
	/**
	 * Push a value at the end of the queue.
	 */
	public void push_back(double x) {
		if (currentLength == size())
			pop_front();
		setElementAt(new Double(x), widx);
		currentLength++;
		sum += x;
		widx++;
		if (widx == size())
			widx = 0;
	}
}
