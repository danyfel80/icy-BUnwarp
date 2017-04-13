# icy-BUnwarp #
B-spline registration adaptation for ICY based on BunwarpJ plugin from ImageJ

## README ##
This README would normally document whatever steps are necessary to get this plugin up and running on Icy.

### What is this repository for? ###
* This repository contains the code of the BUnwarp plugin adapted for working with big images.
* Note: This plugin is an adaptation of the BUnwarpJ plugin available for ImageJ. [Go to BUnwarpJ](http://imagej.net/BUnwarpJ)

### How do I get set up? ###

To set up this plugin in a development environment you should follow these steps 
#### Prerequisites ####
* Have a version of Java JDK superior or equal to 8. [Java download page](http://www.oracle.com/technetwork/java/javase/downloads/index.html).
* Have the latest version of Icy installed on your machine. [Icy download page](http://icy.bioimageanalysis.org/)
* Have Gradle installed to be able to download dependencies, create the eclipse project, and to produce jars as well. [Gradle download page](https://gradle.org/install)
* Have Eclipse IDE installed and configured for working with Icy. [Icy4Eclipse plugin installation manual](http://www.herve.name/pmwiki.php/Main/Icy4Eclipse)

#### Setup ####
1. Make sure to have the environment variable ICY_HOME pointing to the forlder containing icy installation folder.
1. Once downloaded the code, open a terminal on the code and type 
   ```
   cd BUnwarp
   gradle eclipse
   ```
   This code will download and set the project dependencies and create the eclipse project to work with.
1. Open eclipse and import an existing project specifying the location of the project created with gradle.
1. Once the project is open in the workspace you should be able to open icy and find the following plugins.
   * BigImageBUnwarp: Used to register two large images
   * BUnwarpSimple: Adaptation of the original BUnwarpJ plugin from ImageJ
   * LoadBigImage: Loads a large image on icy (either by downsampling or by loading a particular tile)
   * BigImageThresholding: Performs a thresholding filter on a large image.
   * SaveBigImage: Saves a image on icy by tiles (Test plugin)
   * BigImageConvertToTiff: Converts an image file from any format to a large tiff.

### Contributors ###

* Code: Daniel Gonzalez
* Code review: Vannary Meas-Yedid
* Original BUnwarpJ plugin: [Go to BUnwarpJ](http://imagej.net/BUnwarpJ)

### Who do I talk to? ###

* Daniel Felipe Gonzalez Obando: (danyfel80@gmail.com)
