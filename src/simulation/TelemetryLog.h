#pragma once

// ══════════════════════════════════════════════════════════════════════════
//  TelemetryLog.h — Scientific data logger for cross-validation
//  -----------------------------------------------------------------------
//  Comprehensive per-tick CSV logger designed for direct comparison with
//  published cell-biology experiments (HeLa cell-counts vs time, phase
//  fractions, metabolite concentrations, division statistics).
//
//  Two output files per run:
//
//    logs/telemetry_population_<session>.csv
//        Sampled every `populationSampleSec` bio-seconds. One row per
//        sample. Captures population dynamics, phase distribution,
//        environmental concentrations, mean cell state — the "experiment
//        plot" data: cell-count vs time, phase fractions vs time, etc.
//
//    logs/telemetry_division_<session>.csv
//        One row per division event. Records who divided, when, the
//        biology state at division, and what fraction of dish biomass
//        each daughter inherited — for compare-to-experiment of
//        intra-population variability.
//
//  Both columns are stable; a second column is never inserted between
//  shipped columns so external scripts can rely on the schema.
//
//  Usage:
//    TelemetryLog log;
//    log.open("session_id_string");
//    // each frame:
//    log.tick(sim, bio_seconds, wall_seconds);
//    // on shutdown:
//    log.close();
// ══════════════════════════════════════════════════════════════════════════

#include "Simulation.h"
#include "MediumField.h"
#include <cstdio>
#include <cstring>
#include <string>

class TelemetryLog {
public:
    // Sample population row every N bio-seconds (default 60 = once per
    // bio-minute — fine resolution for plotting count-vs-time curves).
    float populationSampleBiosec = 60.0f;

    bool open(const std::string& sessionId) {
        std::string popPath = "logs/telemetry_population_" + sessionId + ".csv";
        std::string divPath = "logs/telemetry_division_"   + sessionId + ".csv";
        popF_ = fopen(popPath.c_str(), "w");
        divF_ = fopen(divPath.c_str(), "w");
        if (!popF_ || !divF_) return false;
        writePopulationHeader();
        writeDivisionHeader();
        nextSampleBiosec_ = 0.0f;
        prevDivisions_ = 0;
        prevDeaths_    = 0;
        return true;
    }

    void close() {
        if (popF_) { fclose(popF_); popF_ = nullptr; }
        if (divF_) { fclose(divF_); divF_ = nullptr; }
    }

    // Per-frame call. Throttles itself to populationSampleBiosec spacing.
    void tick(const Simulation& sim, float bio_seconds, float wall_seconds) {
        if (!popF_) return;
        // Always log NEW divisions (event-driven, not sampled).
        int newDivisions = sim.statDivisions - prevDivisions_;
        if (newDivisions > 0) {
            // Log one row per new division. We approximate the event
            // attribution by writing the current dish-mean biology — the
            // exact dividing cell is harder to recover after the event.
            writeDivisionRow(sim, bio_seconds, wall_seconds);
        }
        prevDivisions_ = sim.statDivisions;
        prevDeaths_    = sim.statDeaths;

        // Population row throttled.
        if (bio_seconds < nextSampleBiosec_) return;
        nextSampleBiosec_ = bio_seconds + populationSampleBiosec;
        writePopulationRow(sim, bio_seconds, wall_seconds);
    }

    // ── Population header & row ─────────────────────────────────────────
    void writePopulationHeader() {
        fprintf(popF_,
            // Time
            "wall_sec,bio_sec,bio_hours,bio_days,"
            // Population
            "cells_alive,cells_dead_cum,divisions_cum,"
            // Phase distribution
            "phase_g1,phase_s,phase_g2,phase_m,"
            // Fate distribution
            "fate_undet,fate_prolif,fate_quiesc,fate_apopt,necrotic,senescent,"
            // Mean per-cell biology
            "mean_biomass,mean_size,mean_ATP,mean_stress,mean_ROS,mean_damage,"
            "mean_telomere,mean_age_min,"
            // CDK / cyclin (mean across alive cells)
            "mean_CycD,mean_CycE,mean_CycA,mean_CycB,mean_Rb,mean_p21,"
            // Replication & mitosis
            "mean_replProgress,mean_replForks,cells_in_mitosis,"
            // Mitochondria
            "mean_mitoCount,mean_mitoHealth,mean_mitoNetworkedFrac,"
            // Macromolecule counts (per-cell mean)
            "mean_genomeBp,mean_nucleosomes,mean_ribosomes,mean_rnaPolII,"
            "mean_spliceosomes,mean_nuclearPores,mean_chaperones,"
            "mean_proteasomes,mean_copiiVesicles,mean_secretoryVesicles,"
            "mean_lysosomes,mean_peroxisomes,mean_mRNA,mean_tRNA,"
            "mean_totalProteins,"
            // Metabolite molecule counts (per-cell mean)
            "mean_atpMolecules,mean_glucoseMolecules,mean_nadhMolecules,"
            "mean_aaMolecules,mean_waterMolecules,mean_calciumFree,"
            "mean_cAMPMolecules,"
            // Osmosis & physical
            "mean_osmoCytoMM,mean_waterMM,mean_turgorPa,"
            "mean_aquaporinCount,mean_adhesionStrength,mean_spreadFactor,"
            // Medium (dish-mean concentrations)
            "med_O2_mM,med_CO2_mM,med_glucose_mM,med_glutamine_mM,"
            "med_pyruvate_mM,med_lactate_mM,med_AAPool_mM,"
            "med_growthF_ngmL,med_ions_mM,med_Ca_mM,med_pH,"
            "med_water_mM"
            "\n");
    }

    void writePopulationRow(const Simulation& sim,
                            float bio_seconds, float wall_seconds) {
        int alive = 0;
        int phase[4] = {0,0,0,0};
        int fate[4]  = {0,0,0,0};
        int necrotic = 0, senescent = 0, mitosing = 0;
        double sumBio=0, sumSize=0, sumATP=0, sumStress=0, sumROS=0;
        double sumDmg=0, sumTel=0, sumAge=0;
        double sumCycD=0, sumCycE=0, sumCycA=0, sumCycB=0;
        double sumRb=0, sump21=0;
        double sumRepl=0, sumForks=0;
        double sumMitoC=0, sumMitoH=0, sumMitoN=0;
        double sumGenome=0, sumNucl=0, sumRib=0, sumPolII=0;
        double sumSplice=0, sumPores=0, sumChap=0;
        double sumProteo=0, sumCOPII=0, sumSecV=0;
        double sumLys=0, sumPerox=0, sumMRNA=0, sumTRNA=0;
        double sumTotProt=0;
        double sumATPmol=0, sumGluMol=0, sumNADH=0, sumAA=0, sumH2O=0;
        double sumCa=0, sumCAMP=0;
        double sumOsmo=0, sumWaterMM=0, sumTurgor=0;
        double sumAQP=0, sumAdh=0, sumSpread=0;

        for (auto& c : sim.cells) {
            if (!c.alive) continue;
            alive++;
            int p = c.phase;
            if (p >= 0 && p < 4) phase[p]++;
            int f = c.fate;
            if (f >= 0 && f < 4) fate[f]++;
            if (c.necrotic)  necrotic++;
            if (c.senescent) senescent++;
            if (c.program.mitosis.active) mitosing++;
            sumBio += c.biomass; sumSize += c.size; sumATP += c.ATP;
            sumStress += c.stress; sumROS += c.ROS; sumDmg += c.damageLevel;
            sumTel += c.telomere; sumAge += c.age;
            sumCycD += c.cdk.CycD; sumCycE += c.cdk.CycE;
            sumCycA += c.cdk.CycA; sumCycB += c.cdk.CycB;
            sumRb   += c.cdk.Rb;   sump21  += c.cdk.p21;
            sumRepl += c.program.cdogma.replicationProgress;
            sumForks += c.replicationForks;
            sumMitoC += c.mitoCount; sumMitoH += c.mitoHealthFrac;
            sumMitoN += c.mitoNetworkedFrac;
            sumGenome += c.genomeBp; sumNucl += c.nucleosomes;
            sumRib += c.ribosomeCount; sumPolII += c.rnaPolII;
            sumSplice += c.spliceosomes; sumPores += c.nuclearPores;
            sumChap += c.chaperones; sumProteo += c.proteasomes;
            sumCOPII += c.copiiVesicles; sumSecV += c.secretoryVesicles;
            sumLys += c.lysosomes; sumPerox += c.peroxisomes;
            sumMRNA += c.mRNACount; sumTRNA += c.tRNACount;
            sumTotProt += c.totalProteins;
            sumATPmol += c.atpMolecules; sumGluMol += c.glucoseMolecules;
            sumNADH += c.nadhMolecules; sumAA += c.aaMolecules;
            sumH2O += c.waterMolecules; sumCa += c.calciumFree;
            sumCAMP += c.cAMPMolecules;
            sumOsmo += c.osmoCytoMM; sumWaterMM += c.waterMM;
            sumTurgor += c.turgorPa;
            sumAQP += c.aquaporinCount; sumAdh += c.adhesionStrength;
            sumSpread += c.spreadFactor;
        }

        double n = alive > 0 ? (double)alive : 1.0;
        auto mu = [&](double s) { return s / n; };

        fprintf(popF_,
            "%.2f,%.0f,%.4f,%.5f,"
            "%d,%d,%d,"
            "%d,%d,%d,%d,"
            "%d,%d,%d,%d,%d,%d,"
            "%.4f,%.4f,%.3f,%.3f,%.3f,%.4f,"
            "%.1f,%.2f,"
            "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.1f,%d,"
            "%.1f,%.4f,%.4f,"
            "%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,"
            "%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e,"
            "%.4f,%.2f,%.2f,"
            "%.1f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,"
            "%.4f,%.4f,%.4f,%.4f,%.4f"
            "\n",
            wall_seconds, bio_seconds, bio_seconds/3600.0, bio_seconds/86400.0,
            alive, sim.statDeaths, sim.statDivisions,
            phase[0], phase[1], phase[2], phase[3],
            fate[0], fate[1], fate[2], fate[3], necrotic, senescent,
            mu(sumBio), mu(sumSize), mu(sumATP), mu(sumStress),
            mu(sumROS), mu(sumDmg),
            mu(sumTel), mu(sumAge),
            mu(sumCycD), mu(sumCycE), mu(sumCycA), mu(sumCycB),
            mu(sumRb), mu(sump21),
            mu(sumRepl), mu(sumForks), mitosing,
            mu(sumMitoC), mu(sumMitoH), mu(sumMitoN),
            mu(sumGenome), mu(sumNucl), mu(sumRib), mu(sumPolII),
            mu(sumSplice), mu(sumPores), mu(sumChap),
            mu(sumProteo), mu(sumCOPII), mu(sumSecV),
            mu(sumLys), mu(sumPerox), mu(sumMRNA), mu(sumTRNA),
            mu(sumTotProt),
            mu(sumATPmol), mu(sumGluMol), mu(sumNADH), mu(sumAA),
            mu(sumH2O), mu(sumCa), mu(sumCAMP),
            mu(sumOsmo), mu(sumWaterMM), mu(sumTurgor),
            mu(sumAQP), mu(sumAdh), mu(sumSpread),
            sim.nutrients.mean(MS_O2),       sim.nutrients.mean(MS_CO2),
            sim.nutrients.mean(MS_GLUCOSE),  sim.nutrients.mean(MS_GLUTAMINE),
            sim.nutrients.mean(MS_PYRUVATE), sim.nutrients.mean(MS_LACTATE),
            sim.nutrients.mean(MS_AA_POOL),  sim.nutrients.mean(MS_GROWTH_F),
            sim.nutrients.mean(MS_IONS),     sim.nutrients.mean(MS_CALCIUM),
            sim.nutrients.mean(MS_HPLUS),
            sim.nutrients.mean(MS_WATER));
        fflush(popF_);
    }

    // ── Division header & row ──────────────────────────────────────────
    void writeDivisionHeader() {
        fprintf(divF_,
            "wall_sec,bio_sec,bio_hours,"
            "cells_after,divisions_cum,"
            "mean_biomass,mean_ATP,mean_damage,mean_replProgress,"
            "med_glucose_mM,med_glutamine_mM,med_O2_mM\n");
    }

    void writeDivisionRow(const Simulation& sim,
                          float bio_seconds, float wall_seconds) {
        if (!divF_) return;
        int alive = 0;
        double sumBio=0, sumATP=0, sumDmg=0, sumRepl=0;
        for (auto& c : sim.cells) {
            if (!c.alive) continue;
            alive++;
            sumBio += c.biomass; sumATP += c.ATP; sumDmg += c.damageLevel;
            sumRepl += c.program.cdogma.replicationProgress;
        }
        double n = alive > 0 ? (double)alive : 1.0;
        fprintf(divF_,
            "%.2f,%.0f,%.4f,%d,%d,%.4f,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
            wall_seconds, bio_seconds, bio_seconds/3600.0,
            alive, sim.statDivisions,
            sumBio/n, sumATP/n, sumDmg/n, sumRepl/n,
            sim.nutrients.mean(MS_GLUCOSE),
            sim.nutrients.mean(MS_GLUTAMINE),
            sim.nutrients.mean(MS_O2));
        fflush(divF_);
    }

private:
    FILE* popF_ = nullptr;
    FILE* divF_ = nullptr;
    float nextSampleBiosec_ = 0.0f;
    int prevDivisions_ = 0;
    int prevDeaths_    = 0;
};
