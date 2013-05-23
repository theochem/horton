// Horton is a development platform for electronic structure methods.
// Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
//
// This file is part of Horton.
//
// Horton is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// Horton is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--


#ifndef HORTON_GRID_RTRANSFORM_H
#define HORTON_GRID_RTRANSFORM_H

class RTransform {
    protected:
        int npoint;

    public:
        RTransform(int npoint);
        virtual ~RTransform() {};
        virtual double radius(double t) = 0;
        virtual double deriv(double t) = 0;
        virtual double inv(double r) = 0;

        void radius_array(double* t, double* r, int n);
        void deriv_array(double* t, double* d, int n);
        void inv_array(double* r, double* t, int n);
        int get_npoint() {return npoint;};
    };


class IdentityRTransform : public RTransform {
    public:
        IdentityRTransform(int npoint): RTransform(npoint) {};
        virtual double radius(double t);
        virtual double deriv(double t);
        virtual double inv(double r);
    };


class LinearRTransform : public RTransform {
    private:
        double rmin, rmax, alpha;
    public:
        LinearRTransform(double rmin, double rmax, int npoint);
        virtual double radius(double t);
        virtual double deriv(double t);
        virtual double inv(double r);

        double get_rmin() {return rmin;};
        double get_rmax() {return rmax;};
        double get_alpha() {return alpha;};
    };


class ExpRTransform : public RTransform {
    private:
        double rmin, rmax, alpha;
    public:
        ExpRTransform(double rmin, double rmax, int npoint);
        virtual double radius(double t);
        virtual double deriv(double t);
        virtual double inv(double r);

        double get_rmin() {return rmin;};
        double get_rmax() {return rmax;};
        double get_alpha() {return alpha;};
    };


class ShiftedExpRTransform : public RTransform {
    private:
        double rmin, rshift, rmax, r0, alpha;
    public:
        ShiftedExpRTransform(double rmin, double rshift, double rmax, int npoint);
        virtual double radius(double t);
        virtual double deriv(double t);
        virtual double inv(double r);

        double get_rmin() {return rmin;};
        double get_rshift() {return rshift;};
        double get_rmax() {return rmax;};
        double get_r0() {return r0;};
        double get_alpha() {return alpha;};
    };


class BakerRTransform : public RTransform {
    private:
        double rmax, scale;
    public:
        BakerRTransform(double rmax, int npoint);
        virtual double radius(double t);
        virtual double deriv(double t);
        virtual double inv(double r);

        double get_rmax() {return rmax;};
        double get_scale() {return scale;};
    };

#endif
