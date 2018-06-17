/*
Minetest
Copyright (C) 2015-2018 paramat
Copyright (C) 2015-2018 kwolekr, Ryan Kwolek <kwolekr@minetest.net>

Fractal formulas from http://www.bugman123.com/Hypercomplex/index.html
by Paul Nylander, and from http://www.fractalforums.com, thank you.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#pragma once

#include "mapgen.h"

class BiomeManager;

extern FlagDesc flagdesc_mapgen_fractal[];

struct MapgenTerrainbrotParams : public MapgenParams
{
	u32 spflags = 0;
	float cave_width = 0.09f;
	s16 large_cave_depth = -33;
	s16 lava_depth = -256;
	s16 dungeon_ymin = -31000;
	s16 dungeon_ymax = 31000;
	u16 iterations = 30;
	u16 rank = 30;
	v3f scale = v3f(4096.0, 512.0, 4096.0);
	v3f offset = v3f(1.52, 0.0, 0.0);

	NoiseParams np_seabed;
	NoiseParams np_filler_depth;
	NoiseParams np_cave1;
	NoiseParams np_cave2;

	NoiseParams np_coord;
	NoiseParams np_polynom;

	float river_min;
	float river_max;
	float river_angle;
	float river_width;
	float river_count;
	float river_bank_coef;
	s16 river_height_min;
	s16 river_height_max;
	u16 river_extra_iterations;
	u16 river_min_iterations;

	MapgenTerrainbrotParams();
	~MapgenTerrainbrotParams() = default;

	void readParams(const Settings *settings);
	void writeParams(Settings *settings) const;
};

class MapgenTerrainbrot : public MapgenBasic
{
public:
	MapgenTerrainbrot(int mapgenid, MapgenTerrainbrotParams *params, EmergeManager *emerge);
	~MapgenTerrainbrot();

	virtual MapgenType getType() const { return MAPGEN_TERRAINBROT; }

	virtual void makeChunk(BlockMakeData *data);
	int getSpawnLevelAtPoint(v2s16 p);
	bool getFractalAtPoint(s16 x, s16 y, s16 z, float* cache, bool& is_water);
	s16 generateTerrain();

private:
	s16 large_cave_depth;
	s16 dungeon_ymin;
	s16 dungeon_ymax;
	//u16 fractal;
	u16 iterations;
	u16 rank;
	v3f scale;
	v3f offset;

	Noise *noise_seabed;
	Noise *noise_coord;
	Noise *noise_polynom;

	float river_min;
	float river_max;
	float river_angle;
	float river_width;
	float river_count;
	float river_bank_coef;
	s16 river_height_min;
	s16 river_height_max;
	u16 river_extra_iterations;
	u16 river_min_iterations;

	void product(float& x, float& y, float a, float b);
	void divide(float& x, float& y, float a, float b);
	void sum(float& x, float& y, float a, float b);

	void polynom(float& x, float& y, s32 seed, u16 my_rank, float xx, float yy, float* cache);
	void rational(float& x, float& y, s32 seed, float xx, float yy, float* cache);
};
