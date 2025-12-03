/**
 * readcross.h
 * 
 * Header file for gas opacity reading and reinterpolation functions
 */

#ifndef READCROSS_H
#define READCROSS_H

// Opacity reinterpolation functions
extern void reinterpolate_all_opacities(void);
extern void cleanup_opacity_cache(void);

// Main opacity reading function
void readcross(char Fname[], double **xsc);
void read_all_opacities(void);

#endif // READCROSS_H

