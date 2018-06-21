/*
 * Copyright 2010-2018 Institut Pasteur.
 * 
 * This file is part of Icy.
 * 
 * Icy is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Icy is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Icy. If not, see <http://www.gnu.org/licenses/>.
 */
package algorithms.danyfel80.command;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;

/**
 * This class is meant to control Registration processes from the command line
 * instead of using Icy's GUI (i. e. headless mode).
 * 
 * @author Daniel Felipe Gonzalez Obando
 */
public class CommandProcessor implements Callable<Void> {
	Map<String, Class<CommandProcess<?>>> commands;

	public CommandProcessor() throws IOException, ClassNotFoundException {
		loadCommandProcessList();
	}

	private void loadCommandProcessList() throws IOException, ClassNotFoundException {
		InputStream in = CommandProcessor.class.getResourceAsStream("/bUnwarpCommandProcessList.txt");
		if (in == null)
			throw new IOException("Could not find process list file");

		commands = new HashMap<>();

		ClassLoader loader = CommandProcessor.class.getClassLoader();
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String line;
		while ((line = br.readLine()) != null) {
			Class<?> loadedClass = loader.loadClass(line);
			if (CommandProcess.class.isAssignableFrom(loadedClass)) {
				@SuppressWarnings("unchecked")
				Class<CommandProcess<?>> loadedProcessClass = (Class<CommandProcess<?>>) loadedClass;
				try {
					commands.put(loadedProcessClass.newInstance().getCommand(), loadedProcessClass);
				} catch (InstantiationException | IllegalAccessException e) {
					throw new ClassNotFoundException("Could not get the command name for class " + loadedProcessClass.getName(),
							e);
				}
			}
		}
	}

	@Override
	public Void call() throws Exception {
		// Open reading stream
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		String line = br.readLine();
		// Read each line
		// Separate parameters by space an convert to lower case
		String[] params = line.split("[ \t]+");
		if (params.length < 1) {
			System.err.println("No command received.");
			return null;
		}
		params = Arrays.stream(params).map(s -> s.toLowerCase()).toArray(String[]::new);

		if (params[0].equals("help")) {
			System.out.println("Commands:");
			commands.entrySet().stream().forEach(k -> {
				try {
					System.out.println(CommandProcessUtil.getCommandDescription(k.getValue()));
				} catch (InstantiationException | IllegalAccessException e) {
					throw new RuntimeException(e);
				}
			});
			return null;
		}

		// Get the process
		Class<CommandProcess<?>> processClass = commands.get(params[0]);
		if (processClass == null) {
			System.err.format("Command %s not recognized.\n", params[0]);
			return null;
		}
		CommandProcess<?> process = processClass.newInstance();

		// Pass parameters
		params = Arrays.stream(params).skip(1).toArray(String[]::new);
		process.setArguments(params);

		// Execute
		Object result = null;
		try {
			result = process.call();
		} catch (Exception e) {
			// When error then print error message
			e.printStackTrace();
			// Print usage when illegal argument is thrown
			if (e instanceof IllegalArgumentException) {
				System.err.println("usage: " + process.toString());
			}
			return null;
		}

		if (result instanceof String)
			System.out.println(result);

		return null;
	}
}
