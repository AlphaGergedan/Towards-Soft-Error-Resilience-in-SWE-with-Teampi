/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * TODO
 */

#ifndef __SWE_SIMPLE_SCENARIOS_H
#define __SWE_SIMPLE_SCENARIOS_H

#include <cmath>

#include "scenarios/Scenario.hpp"

class SWE_HalfDomainDry : public SWE_Scenario
{
public:
    float getBathymetry(float x, float y)
    {
        // if(y > 500 && y < 600)return 20;
        return -250;
    };

    float getWaterHeight(float x, float y)
    {
        // if(y > 500)return 0;
        return (sqrt((x - 250.f) * (x - 250.f) + (y - 250.f) * (y - 250.f)) < 100.f) ? 15.f : 10.0f;
    };

    virtual float endSimulation() { return (float)15; };

    virtual BoundaryType getBoundaryType(Boundary boundary) { return WALL; };

    /** Get the boundary positions
     *
     * @param boundary which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(Boundary boundary)
    {
        if (boundary == BND_LEFT)
            return (float)0;
        else if (boundary == BND_RIGHT)
            return (float)1000;
        else if (boundary == BND_BOTTOM)
            return (float)0;
        else
            return (float)1000;
    };
};

/**
 * Scenario "Radial Dam Break":
 * elevated water in the center of the domain
 */
class SWE_RadialDamBreakScenario : public SWE_Scenario
{
public:
    float getBathymetry(float x, float y) { return 10.f; };

    float getWaterHeight(float x, float y)
    {
        return (sqrt((x - 250.f) * (x - 250.f) + (y - 250.f) * (y - 250.f)) < 100.f) ? 15.f : 10.0f;
    };

    virtual float endSimulation() { return (float)15; };

    virtual BoundaryType getBoundaryType(Boundary boundary) { return OUTFLOW; };

    /** Get the boundary positions
     *
     * @param boundary which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(Boundary boundary)
    {
        if (boundary == BND_LEFT)
            return (float)0;
        else if (boundary == BND_RIGHT)
            return (float)1000;
        else if (boundary == BND_BOTTOM)
            return (float)0;
        else
            return (float)1000;
    };
};

class SWE_RadialBathymetryDamBreakScenario : public SWE_Scenario
{
public:
    float getBathymetry(float x, float y)
    {
        return (std::sqrt((x - 1000.f) * (x - 1000.f) + (y - 1000.f) * (y - 1000.f)) < 300.f)
                   // || (std::sqrt((x - 1000.f) * (x - 1000.f) + (y - 500.f) * (y - 500.f)) < 50.f)
                   ? 20.f
                   : 10.f;
        return -10.f;
    };

    virtual float endSimulation() { return (float)15; };

    virtual BoundaryType getBoundaryType(Boundary edge) { return OUTFLOW; };

    /** Get the boundary positions
     *
     * @param boundary which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(Boundary boundary)
    {
        if (boundary == BND_LEFT)
            return (float)0;
        else if (boundary == BND_RIGHT)
            return (float)2000;
        else if (boundary == BND_BOTTOM)
            return (float)0;
        else
            return (float)2000;
    };

    /**
     * Get the water height at a specific location.
     *
     * @param i_positionX position relative to the origin of the bathymetry grid in x-direction
     * @param i_positionY position relative to the origin of the bathymetry grid in y-direction
     * @return water height (before the initial displacement)
     */
    float getWaterHeight(float x, float y)
    {
        return (sqrt((x - 500.f) * (x - 500.f) + (y - 500.f) * (y - 500.f)) < 300.f) ? 40.f : 15.0f;
    };
};
/**
 * Scenario "Bathymetry Dam Break":
 * uniform water depth, but elevated bathymetry in the centre of the domain
 */
class SWE_BathymetryDamBreakScenario : public SWE_Scenario
{
public:
    float getBathymetry(float x, float y)
    {
        return (std::sqrt((x - 500.f) * (x - 500.f) + (y - 500.f) * (y - 500.f)) < 50.f) ? -255.f : -260.f;
    };

    virtual float endSimulation() { return (float)15; };

    virtual BoundaryType getBoundaryType(Boundary edge) { return OUTFLOW; };

    /** Get the boundary positions
     *
     * @param boundary which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(Boundary boundary)
    {
        if (boundary == BND_LEFT)
            return (float)0;
        else if (boundary == BND_RIGHT)
            return (float)1000;
        else if (boundary == BND_BOTTOM)
            return (float)0;
        else
            return (float)1000;
    };

    /**
     * Get the water height at a specific location.
     *
     * @param i_positionX position relative to the origin of the bathymetry grid in x-direction
     * @param i_positionY position relative to the origin of the bathymetry grid in y-direction
     * @return water height (before the initial displacement)
     */
    float getWaterHeight(float i_positionX, float i_positionY) { return (float)260; }
};

/**
 * Scenario "Sea at Rest":
 * flat water surface ("sea at rest"),
 * but non-uniform bathymetry (id. to "Bathymetry Dam Break")
 * test scenario for "sea at rest"-solution
 */
class SWE_SeaAtRestScenario : public SWE_Scenario
{
public:
    float getWaterHeight(float x, float y)
    {
        return 10.0f;
    };

    float getBathymetry(float x, float y)
    {
        return 0.0f;
    };
    // Added
    virtual float endSimulation() { return (float)15; };
    virtual BoundaryType getBoundaryType(Boundary edge) { return OUTFLOW;};
    float getBoundaryPos(Boundary boundary) {
        if (boundary == BND_LEFT)
            return (float)0;
        else if (boundary == BND_RIGHT)
            return (float)2000;
        else if (boundary == BND_BOTTOM)
            return (float)0;
        else
            return (float)2000;
    }
};

/**
 * Scenario "Splashing Pool":
 * intial water surface has a fixed slope (diagonal to x,y)
 */
class SWE_SplashingPoolScenario : public SWE_Scenario
{
public:
    float getBathymetry(float x, float y) { return -250.f; };

    float getWaterHeight(float x, float y) { return 250.0f + (5.0f - (x + y) / 200); };

    virtual float endSimulation() { return (float)15; };

    virtual BoundaryType getBoundaryType(Boundary edge) { return OUTFLOW; };

    /** Get the boundary positions
     *
     * @param boundary which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(Boundary boundary)
    {
        if (boundary == BND_LEFT)
            return (float)0;
        else if (boundary == BND_RIGHT)
            return (float)2000;
        else if (boundary == BND_BOTTOM)
            return (float)0;
        else
            return (float)2000;
    };
};

/**
 * Scenario "Splashing Cone":
 * bathymetry forms a circular cone
 * intial water surface designed to form "sea at rest"
 * but: elevated water region in the centre (similar to radial dam break)
 */
class SWE_SplashingConeScenario : public SWE_Scenario
{
public:
    float getWaterHeight(float x, float y)
    {
        float r = sqrt((x - 0.5f) * (x - 0.5f) + (y - 0.5f) * (y - 0.5f));
        float h = 4.0f - 4.5f * (r / 0.5f);

        if (r < 0.1f)
            h = h + 1.0f;

        return (h > 0.0f) ? h : 0.0f;
    };

    float getBathymetry(float x, float y)
    {
        float r = sqrt((x - 0.5f) * (x - 0.5f) + (y - 0.5f) * (y - 0.5f));
        return 1.0f + 9.0f * ((r < 0.5f) ? r : 0.5f);
    };

    float waterHeightAtRest() { return 4.0f; };

    float endSimulation() { return 0.5f; };

    virtual BoundaryType getBoundaryType(Boundary boundary) { return OUTFLOW; };
};

#endif
