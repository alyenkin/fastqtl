//FastQTL: Fast and efficient QTL mapper for molecular phenotypes
//Copyright (C) 2015 Olivier DELANEAU, Alfonso BUIL, Emmanouil DERMITZAKIS & Olivier DELANEAU
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "data.h"


void data::runNominal(string fout, double threshold) {

	//0. Prepare genotypes
	vector < double > genotype_sd = vector < double > (genotype_count, 0.0);
	vector < double > phenotype_sd = vector < double > (phenotype_count, 0.0);
	if (covariate_count > 0) {
		LOG.println("\nCorrecting genotypes & phenotypes for covariates");
		covariate_engine->residualize(genotype_orig);
		covariate_engine->residualize(phenotype_orig);
	}
	for (int g = 0 ; g < genotype_count ; g ++) genotype_sd[g] = RunningStat(genotype_orig[g]).StandardDeviation();
	for (int p = 0 ; p < phenotype_count ; p ++) phenotype_sd[p] = RunningStat(phenotype_orig[p]).StandardDeviation();
	normalize(genotype_orig);
	normalize(phenotype_orig);

	//1. Loop over phenotypes
	ofile fdo (fout);
	for (int p = 0 ; p < phenotype_count ; p ++) {

		LOG.println("\nProcessing gene [" + phenotype_id[p] + "]");

		//1.1. Enumerate all genotype-phenotype pairs within cis-window
		vector < int > targetGenotypes, targetDistances;
		for (int g = 0 ; g < genotype_count ; g ++) {
      int cisdistance_tostart;
      int cisdistance_toend;
		  int cisdistance;
		  int distance_startvar_startpheno = genotype_pos[g] - phenotype_start[p];
		  int distance_endvar_startpheno = genotype_end[g] - phenotype_start[p];

      int distance_startvar_endpheno = genotype_pos[g] - phenotype_end[p];
      int distance_endvar_endpheno = genotype_end[g] - phenotype_end[p];

		  // for INVs ignore the span and define the cisdistance
		  // as the distance from the breakpoints to the phenotype_start
		  if (genotype_vartype[g].compare("INV") == 0) {
		    if (abs(distance_startvar_startpheno) <= abs(distance_endvar_startpheno)) {
		      cisdistance_tostart = distance_startvar_startpheno;
        } else {
		      cisdistance_tostart = distance_endvar_startpheno;
        }

        if (abs(distance_startvar_endpheno) <= abs(distance_endvar_endpheno)) {
          cisdistance_toend = distance_startvar_endpheno;
        } else {
          cisdistance_toend = distance_endvar_endpheno;
        }
		  }

		  // for the variants with span (DEL, DUP, MEI), cisdistance_tostart is zero
		  // if the phenotype_start falls within the span, and the distance to
		  // the closest edge otherwise
		  // BNDs get processed here as well, but their END coordinate is the
		  // same as the START coordinate.
		  else {
		    if (distance_startvar_startpheno < 0 && distance_endvar_startpheno > 0) { // if gene is within SV, then cis distance is 0
		      cisdistance_tostart = 0;
		    } else if (distance_startvar_startpheno >= 0) {
		      cisdistance_tostart = distance_startvar_startpheno;
        } else {
		      cisdistance_tostart = distance_endvar_startpheno;
        }

        if (distance_startvar_endpheno < 0 && distance_endvar_endpheno > 0) {
          cisdistance_toend = 0;
        } else if (distance_startvar_endpheno >= 0) {
          cisdistance_toend = distance_startvar_endpheno;
        } else {
          cisdistance_toend = distance_endvar_endpheno;
        }
		  }

      if (cisdistance_tostart > 0 && cisdistance_toend < 0) {
        cisdistance = 0;
      } else if (abs(cisdistance_tostart) < abs(cisdistance_toend)) {
        cisdistance = cisdistance_tostart;
      } else {
        cisdistance = cisdistance_toend;
      }

		  if (abs(cisdistance) <= cis_window) {
		    targetGenotypes.push_back(g);
		    targetDistances.push_back(cisdistance);
		  }
		}
		LOG.println("  * Number of variants in cis = " + sutils::int2str(targetGenotypes.size()));

		//1.2. Nominal pass: scan cis-window & compute statistics
		for (int g = 0 ; g < targetGenotypes.size() ; g ++) {
			double corr = getCorrelation(genotype_orig[targetGenotypes[g]], phenotype_orig[p]);
			double df = sample_count - 2 - covariate_count;
			double tstat2 = getTstat2(corr, df);
			double pval = getPvalueFromTstat2(tstat2, df);
			double slope = getSlope(corr, phenotype_sd[p], genotype_sd[targetGenotypes[g]]);
			double slope_se = abs(slope) / sqrt(tstat2);
			if (pval <= threshold ) {
				fdo << phenotype_id[p];
				fdo << "\t" << genotype_id[targetGenotypes[g]];
				fdo << "\t" << targetDistances[g];
				fdo << "\t" << corr;
        fdo << "\t" << tstat2;
				fdo << "\t" << pval;
				fdo << "\t" << slope;
				fdo << "\t" << slope_se;
				fdo << endl;
			}
		}
		LOG.println("  * Progress = " + sutils::double2str((p+1) * 100.0 / phenotype_count, 1) + "%");
	}
	fdo.close();
}
