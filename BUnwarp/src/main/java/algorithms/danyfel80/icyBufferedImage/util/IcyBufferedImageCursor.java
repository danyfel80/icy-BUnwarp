/*
 * Copyright 2010-2016 Institut Pasteur.
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
package algorithms.danyfel80.icyBufferedImage.util;

import java.util.concurrent.atomic.AtomicBoolean;

import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.TypeUtil;

/**
 * @author Daniel Felipe Gonzalez Obando
 */
public class IcyBufferedImageCursor {

	IcyBufferedImage	plane;
	int[]						planeSize;
	DataType				planeType;

	AtomicBoolean planeChanged;

	Object volumeData;

	/**
	 * 
	 */
	public IcyBufferedImageCursor(IcyBufferedImage vol) {
		this.plane = vol;
		this.planeSize = new int[] { vol.getSizeX(), vol.getSizeY()};
		this.planeType = vol.getDataType_();
		this.volumeData = vol.getDataXYC();
		planeChanged = new AtomicBoolean(false);
	}

	public IcyBufferedImageCursor(Sequence seq, int t, int z) {
		this(seq.getImage(t, z));
	}

	public synchronized double get(int x, int y, int c) {

		switch (planeType) {
		case UBYTE:
		case BYTE:
			return TypeUtil.toDouble(((byte[][]) volumeData)[c][x + y * planeSize[0]], planeType.isSigned());
		case USHORT:
		case SHORT:
			return TypeUtil.toDouble(((short[][]) volumeData)[c][x + y * planeSize[0]], planeType.isSigned());
		case UINT:
		case INT:
			return TypeUtil.toDouble(((int[][]) volumeData)[c][x + y * planeSize[0]], planeType.isSigned());
		case FLOAT:
			return ((float[][]) volumeData)[c][x + y * planeSize[0]];
		case DOUBLE:
			return ((double[][]) volumeData)[c][x + y * planeSize[0]];
		default:
			throw new RuntimeException("Unsupported data type: " + planeType);
		}
	}

	public synchronized void set(int x, int y, int c, double val) {

		switch (planeType) {
		case UBYTE:
		case BYTE:
			((byte[][]) volumeData)[c][x + y * planeSize[0]] = (byte) val;
		break;
		case USHORT:
		case SHORT:
			((short[][]) volumeData)[c][x + y * planeSize[0]] = (short) val;
		break;
		case UINT:
		case INT:
			((int[][]) volumeData)[c][x + y * planeSize[0]] = (int) val;
		break;
		case FLOAT:
			((float[][]) volumeData)[c][x + y * planeSize[0]] = (float) val;
		break;
		case DOUBLE:
			((double[][]) volumeData)[c][x + y * planeSize[0]] = val;
		break;
		default:
			throw new RuntimeException("Unsupported data type");
		}
		planeChanged.set(true);
	}

	public synchronized void setSafe(int x, int y, int c, double val) {

		switch (planeType) {
		case UBYTE:
		case BYTE:
			((byte[][]) volumeData)[c][x + y * planeSize[0]] = (byte) getSafeValue(val);
		break;
		case USHORT:
		case SHORT:
			((short[][]) volumeData)[c][x + y * planeSize[0]] = (short) getSafeValue(val);
		break;
		case UINT:
		case INT:
			((int[][]) volumeData)[c][x + y * planeSize[0]] = (int) getSafeValue(val);
		break;
		case FLOAT:
			((float[][]) volumeData)[c][x + y * planeSize[0]] = (float) getSafeValue(val);
		break;
		case DOUBLE:
			((double[][]) volumeData)[c][x + y * planeSize[0]] = val;
		break;
		default:
			throw new RuntimeException("Unsupported data type");
		}
		planeChanged.set(true);
	}

	private double getSafeValue(double val) {
		return Math.min(Math.max(val, planeType.getMaxValue()), planeType.getMinValue());
	}

	public synchronized void commitChanges() {
		plane.dataChanged();
	}

}
