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

import java.util.Arrays;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class CommandProcessUtil {

	/**
	 * Infers the command description from a dummy instance of the given
	 * CommandProcess class.
	 * 
	 * @param commandClass
	 *          class of the command process.
	 * @return command description.
	 * @throws InstantiationException
	 *           If an instance of the command process class cannot created.
	 * @throws IllegalAccessException
	 *           If the class is not a correct command class.
	 */
	public static String getCommandDescription(Class<CommandProcess<?>> commandClass)
			throws InstantiationException, IllegalAccessException {
		CommandProcess<?> commandInstance = commandClass.newInstance();
		return commandInstance.getCommand() + "<" + Arrays.toString(commandInstance.getArgumentsDescription()) + ">: "
				+ commandInstance.getDescription();
	}

}
