/**
 * 
 */
package algorithms.danyfel80.registration.bunwarp;

/**
 * Enumeration specifying the registration mode
 * @author Daniel Felipe Gonzalez Obando
 */
public enum RegistrationModeEnum {
 FAST(0, "Fast"),
 ACCURATE(1, "Accurate"),
 MONO(2, "Mono");
 
 private final String name;
 private final int number;
 RegistrationModeEnum(int number, String name) {
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
