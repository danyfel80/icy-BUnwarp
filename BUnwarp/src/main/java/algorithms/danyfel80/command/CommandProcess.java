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

import java.util.concurrent.Callable;

/**
 * Represents a command that can be called returning an object of type {@link T}
 * as a result.
 * 
 * @author Daniel Felipe Gonzalez Obando
 *
 * @param <T>
 *          Type of the object returned when the command process is executed
 *          using the {@link #call()} method.
 */
public interface CommandProcess<T> extends Callable<T> {
	
	/**
	 * @return The command recognized to execute this command process.
	 */
	String getCommand();

	/**
	 * @return The name of the command process.
	 */
	String getName();

	/**
	 * @return The description of each argument needed by this command.
	 */
	String[] getArgumentsDescription();

	/**
	 * @return The description of this command process.
	 */
	String getDescription();

	/**
	 * @param args
	 *          Arguments used on the execution of this command process.
	 * @return This command process with the arguments set to {@code args}. This
	 *         allows to chain methods after setting the arguments.
	 */
	CommandProcess<T> setArguments(String[] args);

}
