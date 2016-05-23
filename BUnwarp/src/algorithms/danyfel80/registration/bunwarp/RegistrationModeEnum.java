/**
 * 
 */
package algorithms.danyfel80.registration.bunwarp;

/**
 * Enumeration specifying the registration mode
 * @author Daniel Felipe Gonzalez Obando
 */
public enum RegistrationModeEnum {
 FAST("Fast"),
 ACCURATE("Accurate"),
 MONO("Mono");
 
 private final String name;
 RegistrationModeEnum(String name) {
	 this.name = name;
 }
 
 @Override
 public String toString() {
	 return name;
 }
}
