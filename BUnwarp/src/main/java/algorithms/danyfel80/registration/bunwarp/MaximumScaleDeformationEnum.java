package algorithms.danyfel80.registration.bunwarp;

/**
 * Representation of the minimum scale deformation levels.
 * @author Daniel Felipe Gonzalez Obando
 */
public enum MaximumScaleDeformationEnum {
	VERY_COARSE(0, "Very Coarse"),
	COARSE(1, "Coarse"),
	FINE(2, "Fine"),
	VERY_FINE(3, "Very Fine"),
	SUPER_FINE(4, "Super Fine");
	
	private final int number;
	private final String name;
	private MaximumScaleDeformationEnum(int num, String name) {
		this.number = num;
		this.name = name;
	}
	
	/**
	 * @return the number
	 */
	public int getNumber() {
		return number;
	}
	
	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}
	
	
}