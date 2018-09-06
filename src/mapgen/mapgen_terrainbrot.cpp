/*
Minetest
Copyright (C) 2015-2018 paramat
Copyright (C) 2015-2018 kwolekr, Ryan Kwolek <kwolekr@minetest.net>

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

#include <math.h>

#include "mapgen.h"
#include <cmath>
#include "voxel.h"
#include "noise.h"
#include "mapblock.h"
#include "mapnode.h"
#include "map.h"
#include "content_sao.h"
#include "nodedef.h"
#include "voxelalgorithms.h"
//#include "profiler.h" // For TimeTaker
#include "settings.h" // For g_settings
#include "emerge.h"
#include "dungeongen.h"
#include "cavegen.h"
#include "mg_biome.h"
#include "mg_ore.h"
#include "mg_decoration.h"
#include "mapgen_terrainbrot.h"


FlagDesc flagdesc_mapgen_terrainbrot[] = {
	{NULL,    0}
};

///////////////////////////////////////////////////////////////////////////////////////


MapgenTerrainbrot::MapgenTerrainbrot(int mapgenid, MapgenTerrainbrotParams *params, EmergeManager *emerge)
	: MapgenBasic(mapgenid, params, emerge)
{
	spflags          = params->spflags;
	cave_width       = params->cave_width;
	large_cave_depth = params->large_cave_depth;
	lava_depth       = params->lava_depth;
	dungeon_ymin     = params->dungeon_ymin;
	dungeon_ymax     = params->dungeon_ymax;
	iterations       = params->iterations;
	rank             = params->rank;
	y_scale          = params->y_scale;
	y_offset         = params->y_offset;

	river_min              = params->river_min;
	river_max              = params->river_max;
	river_count            = params->river_count;
	river_angle            = params->river_angle;
	river_width            = params->river_width;
	river_bank_coef        = params->river_bank_coef;
	river_height_min       = params->river_height_min;
	river_height_max       = params->river_height_max;
	river_extra_iterations = params->river_extra_iterations;
	river_min_iterations   = params->river_min_iterations;

	use_cavebrot      = params->use_cavebrot;
	cave_iterations   = params->cave_iterations;
	escape_distance_2 = params->escape_distance * params->escape_distance;

	std::cout << river_min << " " << river_max << " " << river_width << " " << river_count << " " << river_angle << std::endl;

	//// 2D terrain noise
	noise_seabed       = new Noise(&params->np_seabed, seed, csize.X, csize.Z);
	noise_filler_depth = new Noise(&params->np_filler_depth, seed, csize.X, csize.Z);
	noise_coord_x = new Noise(&params->np_coord, seed, csize.X, csize.Z);
	noise_coord_z = new Noise(&params->np_coord, seed+1, csize.X, csize.Z);
	noise_polynom = new Noise*[rank*8];
	noise_cave_coord_re = new Noise(&params->np_cave_coord, seed+4, csize.X, csize.Y + 2, csize.Z);
	noise_cave_coord_im = new Noise(&params->np_cave_coord, seed+7, csize.X, csize.Y + 2, csize.Z);

	for (int i = rank*8 - 1; i >= 0; i--) {
		noise_polynom[i] = new Noise(&params->np_polynom, seed+i, csize.X, csize.Z);
	}

	MapgenBasic::np_cave1 = params->np_cave1;
	MapgenBasic::np_cave2 = params->np_cave2;
}


MapgenTerrainbrot::~MapgenTerrainbrot()
{
	delete noise_seabed;
	delete noise_filler_depth;
	delete noise_coord_x;
	delete noise_coord_z;
	delete noise_cave_coord_re;
	delete noise_cave_coord_im;

	for (int i = rank*8 - 1; i >= 0; i--) {
		delete noise_polynom[i];
	}
	delete noise_polynom;
}


MapgenTerrainbrotParams::MapgenTerrainbrotParams():
	np_seabed       (-1400,9,  v3f(600, 600, 600), 41900, 5, 0.6, 2.0),
	np_filler_depth (0,  1.2,  v3f(150, 150, 150), 261,   3, 0.7, 2.0),
	np_cave1        (0,   12,  v3f(61,  61,  61),  52534, 3, 0.5, 2.0),
	np_cave2        (0,   12,  v3f(67,  67,  67),  10325, 3, 0.5, 2.0),
	np_coord        (0,    1,  v3f(16384, 16384, 16384),  41900, 5, 0.6, 2.0),
	np_polynom      (0,    1,  v3f(12288, 12288, 12288),  41901, 1, 0.6, 2.0),
	np_cave_coord   (0,    1,  v3f(256, 256, 256), 71830, 6, 0.6, 2.0)
{
}


void MapgenTerrainbrotParams::readParams(const Settings *settings)
{
	settings->getFlagStrNoEx("mgterrainbrot_spflags",      spflags, flagdesc_mapgen_fractal);
	settings->getFloatNoEx("mgterrainbrot_cave_width",     cave_width);
	settings->getS16NoEx("mgterrainbrot_large_cave_depth", large_cave_depth);
	settings->getS16NoEx("mgterrainbrot_lava_depth",       lava_depth);
	settings->getS16NoEx("mgterrainbrot_dungeon_ymin",     dungeon_ymin);
	settings->getS16NoEx("mgterrainbrot_dungeon_ymax",     dungeon_ymax);
	settings->getU16NoEx("mgterrainbrot_iterations",       iterations);
	settings->getU16NoEx("mgterrainbrot_rank",             rank);
	settings->getFloatNoEx("mgterrainbrot_y_scale",        y_scale);
	settings->getFloatNoEx("mgterrainbrot_y_offset",       y_offset);

	settings->getNoiseParams("mgterrainbrot_np_filler_depth", np_filler_depth);
	settings->getNoiseParams("mgterrainbrot_np_cave1",        np_cave1);
	settings->getNoiseParams("mgterrainbrot_np_cave2",        np_cave2);
	settings->getNoiseParams("mgterrainbrot_np_coord",        np_coord);
	settings->getNoiseParams("mgterrainbrot_np_polynom",      np_polynom);

	settings->getFloatNoEx("mgterrainbrot_river_min",            river_min);
	settings->getFloatNoEx("mgterrainbrot_river_max",            river_max);
	settings->getFloatNoEx("mgterrainbrot_river_width",          river_width);
	settings->getFloatNoEx("mgterrainbrot_river_count",          river_count);
	settings->getFloatNoEx("mgterrainbrot_river_angle",          river_angle);
	settings->getFloatNoEx("mgterrainbrot_river_bank_coef",      river_bank_coef);
	settings->getS16NoEx("mgterrainbrot_river_height_min",       river_height_min);
	settings->getS16NoEx("mgterrainbrot_river_height_max",       river_height_max);
	settings->getU16NoEx("mgterrainbrot_river_extra_iterations", river_extra_iterations);
	settings->getU16NoEx("mgterrainbrot_river_min_iterations",   river_min_iterations);

	try {
		use_cavebrot = settings->getBool("mgterrainbrot_use_cavebrot");
	} catch (...) {
		//That smells fishy...
	}
	settings->getNoiseParams("mgterrainbrot_np_cavebrot_coord", np_cave_coord);
	settings->getU16NoEx("mgterrainbrot_cavebrot_iterations",   cave_iterations);
	settings->getFloatNoEx("mgterrainbrot_cavebrot_escape",     escape_distance);
}


void MapgenTerrainbrotParams::writeParams(Settings *settings) const
{
	settings->setFlagStr("mgterrainbrot_spflags",      spflags, flagdesc_mapgen_fractal, U32_MAX);
	settings->setFloat("mgterrainbrot_cave_width",     cave_width);
	settings->setS16("mgterrainbrot_large_cave_depth", large_cave_depth);
	settings->setS16("mgterrainbrot_lava_depth",       lava_depth);
	settings->setS16("mgterrainbrot_dungeon_ymin",     dungeon_ymin);
	settings->setS16("mgterrainbrot_dungeon_ymax",     dungeon_ymax);
	settings->setU16("mgterrainbrot_iterations",       iterations);
	settings->setU16("mgterrainbrot_rank",             rank);
	settings->setFloat("mgterrainbrot_y_scale",        y_scale);
	settings->setFloat("mgterrainbrot_y_offset",       y_offset);

	settings->setNoiseParams("mgterrainbrot_np_filler_depth", np_filler_depth);
	settings->setNoiseParams("mgterrainbrot_np_cave1",        np_cave1);
	settings->setNoiseParams("mgterrainbrot_np_cave2",        np_cave2);
	settings->setNoiseParams("mgterrainbrot_np_coord",        np_coord);
	settings->setNoiseParams("mgterrainbrot_np_polynom",      np_polynom);

	settings->setFloat("mgterrainbrot_river_min",            river_min);
	settings->setFloat("mgterrainbrot_river_max",            river_max);
	settings->setFloat("mgterrainbrot_river_width",          river_width);
	settings->setFloat("mgterrainbrot_river_count",          river_count);
	settings->setFloat("mgterrainbrot_river_angle",          river_angle);
	settings->setFloat("mgterrainbrot_river_bank_coef",      river_bank_coef);
	settings->setS16("mgterrainbrot_river_height_min",       river_height_min);
	settings->setS16("mgterrainbrot_river_height_max",       river_height_max);
	settings->setU16("mgterrainbrot_river_extra_iterations", river_extra_iterations);
	settings->setU16("mgterrainbrot_river_min_iterations",   river_min_iterations);

	settings->setBool("mgterrainbrot_use_cavebrot",             use_cavebrot);
	settings->setNoiseParams("mgterrainbrot_np_cavebrot_coord", np_cave_coord);
	settings->setU16("mgterrainbrot_cavebrot_iterations",       cave_iterations);
	settings->setFloat("mgterrainbrot_cavebrot_escape",         escape_distance);
}


/////////////////////////////////////////////////////////////////


int MapgenTerrainbrot::getSpawnLevelAtPoint(v2s16 p)
{
	bool solid_below = false;  // Dry solid node is present below to spawn on
	bool is_water;
	u8 air_count = 0;  // Consecutive air nodes above the dry solid node
	s16 seabed_level = NoisePerlin2D(&noise_seabed->np, p.X, p.Y, seed);
	// Seabed can rise above water_level or might be raised to create dry land
	s16 search_start = MYMAX(seabed_level, water_level + 1);
	if (seabed_level > water_level)
		solid_below = true;

	for (s16 y = search_start; y <= search_start + 128; y++) {
		if (getFractalAtPoint(p.X, y, p.Y, nullptr, 0, is_water)) {  // Fractal node
			solid_below = true;
			air_count = 0;
		} else if (solid_below) {  // Air above solid node
			air_count++;
			// 3 to account for snowblock dust
			if (air_count == 3)
				return y - 2;
		}
	}

	return MAX_MAP_GENERATION_LIMIT;  // Unsuitable spawn point
}


void MapgenTerrainbrot::makeChunk(BlockMakeData *data)
{
	// Pre-conditions
	assert(data->vmanip);
	assert(data->nodedef);
	assert(data->blockpos_requested.X >= data->blockpos_min.X &&
		data->blockpos_requested.Y >= data->blockpos_min.Y &&
		data->blockpos_requested.Z >= data->blockpos_min.Z);
	assert(data->blockpos_requested.X <= data->blockpos_max.X &&
		data->blockpos_requested.Y <= data->blockpos_max.Y &&
		data->blockpos_requested.Z <= data->blockpos_max.Z);

	this->generating = true;
	this->vm   = data->vmanip;
	this->ndef = data->nodedef;
	//TimeTaker t("makeChunk");

	v3s16 blockpos_min = data->blockpos_min;
	v3s16 blockpos_max = data->blockpos_max;
	node_min = blockpos_min * MAP_BLOCKSIZE;
	node_max = (blockpos_max + v3s16(1, 1, 1)) * MAP_BLOCKSIZE - v3s16(1, 1, 1);
	full_node_min = (blockpos_min - 1) * MAP_BLOCKSIZE;
	full_node_max = (blockpos_max + 2) * MAP_BLOCKSIZE - v3s16(1, 1, 1);

	blockseed = getBlockSeed2(full_node_min, seed);

	for (int i = rank*8 - 1; i >= 0; i--) {
		noise_polynom[i]->perlinMap2D(node_min.X, node_min.Z);
	}

	// Generate base terrain, mountains, and ridges with initial heightmaps
	s16 stone_surface_max_y = generateTerrain();

	// Create heightmap
	updateHeightmap(node_min, node_max);

	// Init biome generator, place biome-specific nodes, and build biomemap
	if (flags & MG_BIOMES) {
		biomegen->calcBiomeNoise(node_min);
		generateBiomes();
	}

	if (use_cavebrot) {
		cavebrot();
	} else {
		if (flags & MG_CAVES) {
			// Generate tunnels
			generateCavesNoiseIntersection(stone_surface_max_y);
			// Generate large randomwalk caves
			generateCavesRandomWalk(stone_surface_max_y, large_cave_depth);
		}
	}

	if ((flags & MG_DUNGEONS) && full_node_min.Y >= dungeon_ymin &&
			full_node_max.Y <= dungeon_ymax)
		generateDungeons(stone_surface_max_y);

	// Generate the registered decorations
	if (flags & MG_DECORATIONS)
		m_emerge->decomgr->placeAllDecos(this, blockseed, node_min, node_max);

	// Generate the registered ores
	m_emerge->oremgr->placeAllOres(this, blockseed, node_min, node_max);

	// Sprinkle some dust on top after everything else was generated
	if (flags & MG_BIOMES)
		dustTopNodes();

	//printf("makeChunk: %dms\n", t.stop());

	updateLiquid(&data->transforming_liquid, full_node_min, full_node_max);

	if (flags & MG_LIGHT)
		calcLighting(node_min - v3s16(0, 1, 0), node_max + v3s16(0, 1, 0),
			full_node_min, full_node_max);

	//setLighting(node_min - v3s16(1, 0, 1) * MAP_BLOCKSIZE,
	//			node_max + v3s16(1, 0, 1) * MAP_BLOCKSIZE, 0xFF);

	this->generating = false;
}


void MapgenTerrainbrot::product(float& x, float& y, float a, float b) {
	float x2 = x*a - y*b;
	y = x*b + y*a;
	x = x2;
}

void MapgenTerrainbrot::sum(float& x, float& y, float a, float b) {
	x += a;
	y += b;
}

void MapgenTerrainbrot::divide(float& x, float& y, float a, float b) {
	float p = a*a + b*b;
	product(x, y, a, -b);
	x /= p;
	y /= p;
}

void MapgenTerrainbrot::polynom(float& x, float& y, s32 my_seed, u16 my_rank, float xx, float yy, Noise** cache, u32 index2d) {
	float cx = 1, cy = 0;
	float bx = x, by = y;
	x = 0;
	y = 0;
	for (u16 i = 0; i <= my_rank; ++i) {
		float ax, ay;
		if (cache) {
			ax = cache[2*i]->result[index2d] * (i+1) / 4;
			ay = cache[2*i+1]->result[index2d] * (i+1) / 4;
		} else {
			ax = NoisePerlin2D(&noise_polynom[2*i]->np, xx, yy, my_seed + 2*i) * (i+1) / 4;
			ay = NoisePerlin2D(&noise_polynom[2*i+1]->np, xx, yy, my_seed + 2*i+1) * (i+1) / 4;
		}
		product(ax, ay, cx, cy);
		sum(x, y, ax, ay);
		product(cx, cy, bx, by);
	}
}

void MapgenTerrainbrot::rational(float& x, float& y, s32 my_seed, float xx, float yy, Noise** cache, u32 index2d) {
	float x2 = x;
	float y2 = y;
	if (cache) {
		polynom(x, y, my_seed, rank, xx, yy, cache, index2d);
		polynom(x2, y2, my_seed + rank*2 + 2, rank-2, xx, yy, &(cache[rank*2 + 2]), index2d);
	} else {
		polynom(x, y, my_seed, rank, xx, yy, nullptr, 0);
		polynom(x2, y2, my_seed + rank*2 + 2, rank-2, xx, yy, nullptr, 0);
	}
	divide(x, y, x2, y2);
}

void MapgenTerrainbrot::cave_rational(float& x, float& y, s32 my_seed, float xx, float yy, Noise** cache, u32 index2d) {
	float x2 = x;
	float y2 = y;
	if (cache) {
		polynom(x, y, my_seed, rank-2, xx, yy, cache, index2d);
		polynom(x2, y2, my_seed + rank*2 + 2, rank, xx, yy, &(cache[rank*2 + 2]), index2d);
	} else {
		polynom(x, y, my_seed, rank-2, xx, yy, nullptr, 0);
		polynom(x2, y2, my_seed + rank*2 + 2, rank, xx, yy, nullptr, 0);
	}
	divide(x, y, x2, y2);
}

bool MapgenTerrainbrot::getFractalAtPoint(s16 x, s16 y, s16 z, Noise** cache, u32 index2d, bool& is_water)
{
	float xx = x, zz = z, cx, cy, cz, ox, oz;
	is_water = false;

	if (cache) {
		cx = noise_coord_x->result[index2d];
		cz = noise_coord_z->result[index2d];
	} else {
		cx = NoisePerlin2D(&noise_coord_x->np, x, z, seed);
		cz = NoisePerlin2D(&noise_coord_z->np, x, z, seed + 1);
	}

	cy = (float)y / y_scale - y_offset;

	float nx = 0.0f;
	float nz = 0.0f;

	float aax, aay;
	aax = cx;
	aay = cz;
	if (cache) {
		rational(aax, aay, seed + rank*4, xx, zz, &(cache[rank*4]), index2d);
	} else {
		rational(aax, aay, seed + rank*4, xx, zz, nullptr, 0);
	}
	ox = cx;
	oz = cz;

	bool count_river = y <= river_height_max && y >= river_height_min;

	u16 iter;
	for (iter = 0; iter < (count_river ? iterations + river_extra_iterations : iterations); iter++) {
		nx = ox;
		nz = oz;
		rational(ox, oz, seed, xx, zz, cache, index2d);
		sum(nx, nz, aax, aay);
		float dist = exp(sqrt((nx-ox)*(nx-ox) + (nz-oz)*(nz-oz)) * cy) * 0.7;
		nx *= dist;
		nz *= dist;

		if (nx * nx + nz * nz > 4.0f || nx != nx || nz != nz) {
			if (iter < iterations) {
				return false;
			} else {
				break;
			}
		}

		ox = nx;
		oz = nz;
	}
	iter -= iterations;

	if (count_river) {
		float mag = sqrt(nx*nx + nz*nz);
		if (mag > river_min) {
			float angle = fabs(fmodf((atan2(nz, nx) + M_PI)*river_count, M_PI*2) - (river_angle + river_width));
			if (angle - (mag - river_min)*river_bank_coef < river_width) {
				if (mag < river_max && iter >= river_min_iterations) {
					is_water = true;
				} else {
					return false;
				}
			}
		}
	}

	return true;
}


s16 MapgenTerrainbrot::generateTerrain()
{
	MapNode n_air(CONTENT_AIR);
	MapNode n_stone(c_stone);
	MapNode n_water(c_water_source);
	MapNode n_river_water(c_river_water_source);

	s16 stone_surface_max_y = -MAX_MAP_GENERATION_LIMIT;
	u32 index2d = 0;
	bool is_water;
	u16 x_size = node_max.X - node_min.X + 1;
	u16 y_size = node_max.Y - node_min.Y + 3;

	noise_seabed->perlinMap2D(node_min.X, node_min.Z);
	noise_coord_x->perlinMap2D(node_min.X, node_min.Z);
	noise_coord_z->perlinMap2D(node_min.X, node_min.Z);


	float *perlin_cache = new float[rank*8 * x_size];
	bool waters_x = false;
	bool *waters_y = new bool[x_size];
	bool *waters_z = new bool[x_size * y_size];

	for (s16 z = node_min.Z; z <= node_max.Z; z++) {

		for (s16 y = node_min.Y - 1; y <= node_max.Y + 1; y++) {
			u32 vi = vm->m_area.index(node_min.X, y, z);
			s16 y_ = y - (node_min.Y - 1);
			for (s16 x = node_min.X; x <= node_max.X; x++, vi++, index2d++) {
				s16 x_ = x - node_min.X;
				if (vm->m_data[vi].getContent() == CONTENT_IGNORE) {
					s16 seabed_height = noise_seabed->result[index2d];

					if (y <= seabed_height || getFractalAtPoint(x, y, z,
						noise_polynom, index2d, is_water)) {
						if (is_water && (waters_y[x_] || waters_x || waters_z[y_*x_size + x_])) {
							vm->m_data[vi] = n_river_water;
						} else {
							vm->m_data[vi] = n_stone;
						}
						if (y > stone_surface_max_y)
							stone_surface_max_y = y;
					} else if (y <= water_level) {
						vm->m_data[vi] = n_water;
					} else {
						vm->m_data[vi] = n_air;
					}
					waters_x = is_water;
					waters_y[x_] = is_water;
					waters_z[y_*x_size + x_] = is_water;
				}
			}
			index2d -= ystride;
		}
		index2d += ystride;
	}
	delete perlin_cache;
	delete waters_y;
	delete waters_z;

	return stone_surface_max_y;
}

bool MapgenTerrainbrot::getCavebrotAtPoint(float re, float im, u32 index2d)
{
	float rr, ii;
	rr = re;
	ii = im;
	for (u16 iter = 0; iter < cave_iterations; iter++) {
		cave_rational(rr, ii, seed, 0, 0, noise_polynom, index2d);
		sum(rr, ii, re, im);
		if (rr*rr + ii*ii > escape_distance_2 || rr != rr || ii != ii) {
			return false;
		}
	}
	return true;	
}

void MapgenTerrainbrot::cavebrot()
{
	MapNode n_air(CONTENT_AIR);
	MapNode n_stone(c_stone);
	MapNode n_water(c_water_source);

	u32 index2d = 0;
	u32 index3d = 0;

	noise_cave_coord_re->perlinMap3D(node_min.X, node_min.Y - 1, node_min.Z);
	noise_cave_coord_im->perlinMap3D(node_min.X, node_min.Y - 1, node_min.Z);

	for (s16 z=node_min.Z; z<=node_max.Z; z++) {
		for (s16 y=node_min.Y - 1; y<=node_max.Y + 1; y++) {
			u32 vi = vm->m_area.index(node_min.X, y, z);
			for (s16 x=node_min.X; x<=node_max.X; x++, vi++, index3d++, index2d++) {

				content_t c = vm->m_data[vi].getContent();
				if (ndef->get(c).is_ground_content) {
					if (!getCavebrotAtPoint(noise_cave_coord_re->result[index3d],
					noise_cave_coord_im->result[index3d], index2d)) {
						vm->m_data[vi] = n_air;
						vm->m_flags[vi] |= VMANIP_FLAG_CAVE;
					}
				}
			}
			index2d -= ystride;
		}
		index2d += ystride;
	}
}
