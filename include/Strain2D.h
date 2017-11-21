/*
* File:   Strain2D.h
* Author: justin
*
* Created on June 7, 2015, 10:14 PM
*/

#ifndef STRAIN2D_H
#define	STRAIN2D_H

#include "Array2D.h"
#include "Image2D.h"
#include "ROI2D.h"
#include "Data2D.h"
#include <math.h>

namespace ncorr {

	namespace details {
		class Strain2D_nlinfo_interpolator;
	}

	class Strain2D final {
		// ---------------------------------------------------------------------------//
		// Strain2D is a class for 2D strains. ---------------------------------------//
		// ---------------------------------------------------------------------------//
	public:
		typedef Data2D::difference_type                             difference_type;
		typedef Data2D::coords                                               coords;
		typedef details::Strain2D_nlinfo_interpolator           nlinfo_interpolator;

		// Rule of 5 and destructor ----------------------------------------------//    
		Strain2D() noexcept = default;
		Strain2D(const Strain2D&) = default;
		Strain2D(Strain2D&&) noexcept = default;
		Strain2D& operator=(const Strain2D&) = default;
		Strain2D& operator=(Strain2D&&) = default;
		~Strain2D() noexcept = default;

		// Additional constructors -----------------------------------------------//
		Strain2D(Array2D<double> eyy, Array2D<double> exy, Array2D<double> exx, const ROI2D &roi, difference_type scalefactor) : // r-value
			eyy(std::move(eyy), roi, scalefactor), exy(std::move(exy), roi, scalefactor), exx(std::move(exx), roi, scalefactor) {
			create_e1_e2();
		}


		// Static factory methods ------------------------------------------------//
		static Strain2D load(std::ifstream&);

		// Operators interface ---------------------------------------------------//
		friend std::ostream& operator<<(std::ostream&, const Strain2D&);
		friend void imshow(const Strain2D& strain, difference_type delay = -1) {
			// Just show each separately for now. If you combine into one buffer, you
			// must scale each as their ranges might be different.
			imshow(strain.eyy, delay);
			imshow(strain.exy, delay);
			imshow(strain.exx, delay);
		}

		friend bool isequal(const Strain2D&, const Strain2D&);
		friend void save(const Strain2D&, std::ofstream&);

		// Access ----------------------------------------------------------------//
		// Note that Strain2D is immutable, so all access is const.    
		difference_type data_height() const { return eyy.data_height(); }
		difference_type data_width() const { return eyy.data_width(); }
		const Data2D& get_eyy() const { return eyy; }
		const Data2D& get_exy() const { return exy; }
		const Data2D& get_exx() const { return exx; }
		const Data2D& get_e1() const { return e1; }
		const Data2D& get_e2() const { return e2; }
		const ROI2D& get_roi() const { return eyy.get_roi(); }
		difference_type get_scalefactor() const { return eyy.get_scalefactor(); }

		// Interpolator ----------------------------------------------------------//
		nlinfo_interpolator get_nlinfo_interpolator(difference_type, INTERP) const;

		// Utility ---------------------------------------------------------------//
		std::string size_string() const { return eyy.size_string(); }
		std::string size_2D_string() const { return eyy.size_2D_string(); }

	private:
		Data2D eyy;   // Immutable - Data2D has pointer semantics
		Data2D exy;   // Immutable - Data2D has pointer semantics
		Data2D exx;   // Immutable - Data2D has pointer semantics
		Data2D e1;
		Data2D e2;
		void create_e1_e2() {

			Array2D<double> e1_array = *new Array2D<double>(eyy.get_array());
			Array2D<double> e2_array = *new Array2D<double>(eyy.get_array());

			auto i_eyy = eyy.get_array().begin();
			auto i_exx = exx.get_array().begin();
			auto i_exy = exy.get_array().begin();
			auto i_e1 = e1_array.begin();
			auto i_e2 = e2_array.begin();

			while (i_eyy != eyy.get_array().end()) {

				double d_eyy = *i_eyy;
				double d_exy = *i_exy;
				double d_exx = *i_exx;
				double det = d_eyy*d_exx - d_exy*d_exy;
				double tr = d_exx + d_eyy;
				double eigenfirst = tr / 2 + std::sqrt(std::pow(d_eyy-d_exx, 2.0) / 4 + std::pow(d_exy, 2.0));
				double eigensecond = tr / 2 - std::sqrt(std::pow(tr, 2.0) / 4 + std::pow(d_exy, 2.0));
				*i_e1 = std::abs(eigenfirst)>std::abs(eigensecond)? eigenfirst: eigensecond;
				*i_e2 = std::abs(eigenfirst)<std::abs(eigensecond)? eigenfirst: eigensecond;
				++i_eyy;
				++i_exx;
				++i_exy;
				++i_e1;
				++i_e2;

			}

			e1 = *new Data2D(e1_array, eyy.get_roi(), eyy.get_scalefactor());
			e2 = *new Data2D(e2_array, eyy.get_roi(), eyy.get_scalefactor());

		}
	};

	namespace details {
		class Strain2D_nlinfo_interpolator final {
		public:
			typedef Strain2D::difference_type                   difference_type;
			typedef Strain2D::coords                                     coords;

			friend Strain2D;

			// Rule of 5 and destructor --------------------------------------//
			Strain2D_nlinfo_interpolator() noexcept = default;
			Strain2D_nlinfo_interpolator(const Strain2D_nlinfo_interpolator&) = default;
			Strain2D_nlinfo_interpolator(Strain2D_nlinfo_interpolator&&) = default;
			Strain2D_nlinfo_interpolator& operator=(const Strain2D_nlinfo_interpolator&) = default;
			Strain2D_nlinfo_interpolator& operator=(Strain2D_nlinfo_interpolator&&) = default;
			~Strain2D_nlinfo_interpolator() noexcept = default;

			// Additional Constructors ---------------------------------------//            
			Strain2D_nlinfo_interpolator(const Strain2D &strain, difference_type region_idx, INTERP interp_type) :
				eyy_interp(strain.get_eyy().get_nlinfo_interpolator(region_idx, interp_type)),
				exy_interp(strain.get_exy().get_nlinfo_interpolator(region_idx, interp_type)),
				exx_interp(strain.get_exx().get_nlinfo_interpolator(region_idx, interp_type)) { }

			// Access methods ------------------------------------------------//
			std::tuple<double, double, double> operator()(double p1, double p2) const { return std::make_tuple(eyy_interp(p1, p2), exy_interp(p1, p2), exx_interp(p1, p2)); }
			std::tuple<const Array2D<double>&, const Array2D<double>&, const Array2D<double>&> first_order(double p1, double p2) const {
				std::cout << "changed method called";
				return std::tuple<const Array2D<double>&, const Array2D<double>&, const Array2D<double>&>(eyy_interp.first_order(p1, p2), exy_interp.first_order(p1, p2), exx_interp.first_order(p1, p2));
			}

		private:
			Data2D::nlinfo_interpolator eyy_interp; // must have copy of interpolator
			Data2D::nlinfo_interpolator exy_interp; // must have copy of interpolator
			Data2D::nlinfo_interpolator exx_interp; // must have copy of interpolator
		};
	}

}

#endif	/* STRAIN2D_H */

