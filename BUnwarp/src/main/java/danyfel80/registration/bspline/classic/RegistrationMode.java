package danyfel80.registration.bspline.classic;

public enum RegistrationMode {
	FAST(0, "Fast"),
	ACCURATE(1, "Accurate"),
	MONO(2, "Mono");

	private final String name;
	private final int number;

	RegistrationMode(int number, String name) {
		this.name = name;
		this.number = number;
	}

	public int getNumber() {
		return number;
	}

	@Override
	public String toString() {
		return name;
	}
}
