#ifndef GRAINITACALOSDACTION_C
#define GRAINITACALOSDACTION_C 1
#endif

/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "GrainitaCaloSDAction.h"
//#include "detectorSegmentations/FCCSWModularGridRhoPhiTheta_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "DD4hep/Segmentations.h"
#include "DDG4/Factories.h"
#include "DDG4/Geant4GeneratorAction.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4SensDetAction.inl"

#include "G4EmProcessSubType.hh"
#include "G4OpticalPhoton.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHandle.hh"
#include "G4VProcess.hh"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>


// #define DEBUG

namespace dd4hep {
namespace sim {

  template <>
  Geant4SensitiveAction<GrainitaCaloSDData>::Geant4SensitiveAction(Geant4Context* ctxt, const std::string& nam,
                                                                    DetElement det, Detector& desc)
      : Geant4Sensitive(ctxt, nam, det, desc), m_collectionName(), m_collectionID(0) {
    declareProperty("ReadoutName", m_readoutName);
    declareProperty("CollectionName", m_collectionName);
    declareProperty("RawCollectionName", m_userData.rawCollectionName);
    declareProperty("useLightResponseFunction", m_userData.useLightResponseFunction);
    declareProperty("responseFuncSlope", m_userData.slope);
    declareProperty("responseFuncIntersect", m_userData.intersect);
    declareProperty("responseFuncX0", m_userData.x0);
    declareProperty("responseFuncAttLength", m_userData.AttLength);
    declareProperty("responseFuncNorm", m_userData.norm);
    declareProperty("neighborCellSize", m_userData.neighborCellSize);
    declareProperty("fiberAttenuationLength", m_userData.fiberAttenuationLength);
    declareProperty("outerRadius", m_userData.outerRadius);


    InstanceCount::increment(this);

    m_hitCreationMode = HitCreationFlags::DETAILED_MODE;
  }

  // Function template specialization of Geant4SensitiveAction class.
  // Define actions
  template <>
  void Geant4SensitiveAction<GrainitaCaloSDData>::initialize() {
    m_userData.sensitive = this;
    m_hitCreationMode = HitCreationFlags::DETAILED_MODE;
  }

  // Function template specialization of Geant4SensitiveAction class.
  // Define collections created by this sensitivie action object
  template <>
  void Geant4SensitiveAction<GrainitaCaloSDData>::defineCollections() {
    m_collectionID = defineCollection<Geant4Calorimeter::Hit>(m_collectionName);
    m_userData.rawCollectionID = defineCollection<Geant4Calorimeter::Hit>(m_userData.rawCollectionName);
  }

  // Function template specialization of Geant4SensitiveAction class.
  // Method that accesses the G4Step object at each track step.
  template <>
  bool Geant4SensitiveAction<GrainitaCaloSDData>::process(const G4Step* aStep, G4TouchableHistory* /*history*/) {

#ifdef DEBUG
    std::cout << "-------------------------------" << std::endl;
    std::cout << "--> GrainitaCalo: track info: " << std::endl;
    std::cout << "----> Track #: " << aStep->GetTrack()->GetTrackID() << " "
              << "Step #: " << aStep->GetTrack()->GetCurrentStepNumber() << " "
              << "Volume: " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << " "
              << std::endl;
    std::cout << "--> GrainitaCalo:: position info(mm): " << std::endl;
    std::cout << "----> x: " << aStep->GetPreStepPoint()->GetPosition().x()
              << " y: " << aStep->GetPreStepPoint()->GetPosition().y()
              << " z: " << aStep->GetPreStepPoint()->GetPosition().z() << std::endl;
    std::cout << "--> GrainitaCalo: particle info: " << std::endl;
    std::cout << "----> Particle " << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << " "
              << "Dep(MeV) " << aStep->GetTotalEnergyDeposit() << " "
              << "Mat " << aStep->GetPreStepPoint()->GetMaterial()->GetName() << " "
              << "Vol " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << " " << std::endl;
#endif

    G4Track* track = aStep->GetTrack();
    if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
      const G4VProcess* creator = track->GetCreatorProcess();
      const bool isCerenkov = creator && creator->GetProcessSubType() == G4EmProcessSubType::fCerenkov;

#ifdef DEBUG
      std::cout << "--> GrainitaCalo: optical photon from "
                << (creator ? creator->GetProcessName() : "unknown")
                << " isCerenkov=" << isCerenkov << std::endl;
#endif

      track->SetTrackStatus(fStopAndKill);
      return true;
    }

    // auto decoder = m_sensitive.readout().idSpec().decoder();
    auto decoder = m_segmentation->decoder();
    auto VolID = volumeID(aStep);
    m_userData.norm = 1. / std::exp(-1. * m_userData.x0 / m_userData.AttLength);


#ifdef DEBUG
    auto SystemID = decoder->get(VolID, "system");
    auto StaveID = decoder->get(VolID, "stave");
    auto SectorID = decoder->get(VolID, "sector");
    std::cout<< "--> Volume ID: "<<SystemID<<"  "<<StaveID<<"  "<<SectorID<<std::endl;
#endif

    // G4TouchableHandle theTouchable = aStep->GetPreStepPoint()->GetTouchableHandle();
    G4ThreeVector global = (aStep->GetPreStepPoint()->GetPosition() + aStep->GetPostStepPoint()->GetPosition() )/2.;
    dd4hep::Position glob(global.x() * dd4hep::millimeter / CLHEP::millimeter,
                          global.y() * dd4hep::millimeter / CLHEP::millimeter,
                          global.z() * dd4hep::millimeter / CLHEP::millimeter);

    auto cellID = m_segmentation->cellID(glob, glob, VolID);
    auto hitpos_dd4hep = m_segmentation->position(cellID); // in cm
    G4ThreeVector HitCellPos(hitpos_dd4hep.x()*dd4hep::centimeter/dd4hep::millimeter, hitpos_dd4hep.y()*dd4hep::centimeter/dd4hep::millimeter, hitpos_dd4hep.z()*dd4hep::centimeter/dd4hep::millimeter );

    Geant4HitCollection* rawColl = collection(m_userData.rawCollectionID);
    const G4double rawStepE = aStep->GetTotalEnergyDeposit();
    Geant4Calorimeter::Hit* rawHit = rawColl->findByKey<Geant4Calorimeter::Hit>(cellID);
    if (!rawHit) {
      rawHit = new Geant4Calorimeter::Hit();
      rawHit->cellID = cellID;
      rawHit->position = HitCellPos;
      rawHit->energyDeposit = rawStepE;
      rawColl->add(cellID, rawHit);
    } else {
      rawHit->energyDeposit += rawStepE;
    }

    Geant4Calorimeter::Hit::Contribution rawContrib;
    rawContrib.trackID = aStep->GetTrack()->GetTrackID();
    rawContrib.pdgID = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
    rawContrib.deposit = rawStepE;
    rawContrib.time = aStep->GetPreStepPoint()->GetGlobalTime();
    rawContrib.x = global.x();
    rawContrib.y = global.y();
    rawContrib.z = global.z();
    rawHit->truth.emplace_back(rawContrib);

    // auto modularSeg = dynamic_cast<const dd4hep::DDSegmentation::FCCSWModularGridRhoPhiTheta_k4geo*>(m_segmentation->segmentation);
    auto phiThetaSeg = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo*>(m_segmentation->segmentation);
    const int phiIndex = decoder->index("phi");
    const int thetaIndex = decoder->index("theta");
    const int phiBins = phiThetaSeg ? phiThetaSeg->phiBins() : 0;
    // cellID() encodes atan2(phi) in the periodic interval whose first bin contains -pi.
    const int firstPhiID = phiThetaSeg
                               ? static_cast<int>(std::floor(
                                     (-M_PI + 0.5 * phiThetaSeg->gridSizePhi() - phiThetaSeg->offsetPhi()) /
                                     phiThetaSeg->gridSizePhi()))
                               : 0;
    const int currentPhiID = static_cast<int>(decoder->get(cellID, phiIndex));
    const int currentThetaID = static_cast<int>(decoder->get(cellID, thetaIndex));
    const int neighborSize = std::max(1, m_userData.neighborCellSize);
    const int neighborRadius = neighborSize / 2;

    std::vector<CellID> cellIDvec;
    std::vector<G4ThreeVector> cellPosVec;
    std::vector<G4double> responseVec;

    auto dd4hepPositionToG4 = [](const dd4hep::DDSegmentation::Vector3D& pos) {
      return G4ThreeVector(pos.X * dd4hep::centimeter / dd4hep::millimeter,
                           pos.Y * dd4hep::centimeter / dd4hep::millimeter,
                           pos.Z * dd4hep::centimeter / dd4hep::millimeter);
    };

    auto addCell = [&](CellID id) {
      auto pos = m_segmentation->position(id);
      cellIDvec.push_back(id);
      cellPosVec.push_back(dd4hepPositionToG4(pos));
    };

    addCell(cellID);
    if(m_userData.useLightResponseFunction){
      for (int dPhi = -neighborRadius; dPhi <= neighborRadius; ++dPhi) {
        for (int dTheta = -neighborRadius; dTheta <= neighborRadius; ++dTheta) {
          if (dPhi == 0 && dTheta == 0) {
            continue;
          }

          int neighborPhiID = currentPhiID + dPhi;
          if (phiBins > 0) {
            const int relativePhiID = neighborPhiID - firstPhiID;
            neighborPhiID = firstPhiID + (relativePhiID % phiBins + phiBins) % phiBins;
          } else {
            std::cout << "Error: phiBins is not well defined: " << phiBins << ". Cannot apply periodic boundary conditions." << std::endl;
            continue;
          }

          const int neighborThetaID = currentThetaID + dTheta;

          CellID neighborCellID = cellID;
          decoder->set(neighborCellID, phiIndex, neighborPhiID);
          decoder->set(neighborCellID, thetaIndex, neighborThetaID);
          addCell(neighborCellID);
        }
      }
    }

    auto transverseDistance = [&](CellID id, const G4ThreeVector& cellPos) {
      G4ThreeVector axis = cellPos.unit();
      // if (modularSeg) {
      //   auto fiberDir = modularSeg->fiberDirection(id);
      //   axis = G4ThreeVector(fiberDir.X, fiberDir.Y, fiberDir.Z).unit();
      // }
      G4ThreeVector rel = global - cellPos;
      G4double transverse2 = rel.mag2() - rel.dot(axis) * rel.dot(axis);

      // std::cout<<""<<"CellID "<<id<<" cellPos ("<<cellPos.x()<<", "<<cellPos.y()<<", "<<cellPos.z()<<") "
      //     <<" global ("<<global.x()<<", "<<global.y()<<", "<<global.z()<<") "
      //     <<" rel ("<<rel.x()<<", "<<rel.y()<<", "<<rel.z()<<") "
      //     <<" axis ("<<axis.x()<<", "<<axis.y()<<", "<<axis.z()<<") "
      //     <<" transverse distance: "<<std::sqrt(std::max(0., transverse2))<<" mm"<<std::endl;

      return std::sqrt(std::max(0., transverse2));
    };

    responseVec.reserve(cellIDvec.size());
    for (std::size_t i = 0; i < cellIDvec.size(); ++i) {
      if(m_userData.useLightResponseFunction){
        double distance = transverseDistance(cellIDvec[i], cellPosVec[i]);
        double response = m_userData.lightResponse(distance);
        // std::cout<<"  Response for cell "<<cellIDvec[i]<<": distance: "<<distance<<", response: "<<response<<std::endl;
        // std::cout<<std::endl;
        responseVec.push_back(response);
      }
      else responseVec.push_back(1.);
    }

    // //Normalize the response in fibers
    // double responseSum = std::accumulate(responseVec.begin(), responseVec.end(), 0.0);
    // std::vector<G4double> responseVec_Norm;
    // for (auto &r : responseVec) {
    //  double r_norm = r / responseSum;
    //  responseVec_Norm.push_back(r_norm);
    // }

#ifdef DEBUG
    auto phiID = m_segmentation->decoder()->get(cellID, "phi");
    auto thetaID = m_segmentation->decoder()->get(cellID, "theta");
    auto rhoID = m_segmentation->decoder()->get(cellID, "rho");
    std::cout<<"--> Step global position: ("<<global.x()<<", "<<global.y()<<", "<<global.z()<<") ";
    std::cout<<" (theta, phi, rho) = "<<"("<<global.theta()<<", "<<global.phi()<<", "<<global.mag()<<") "<<std::endl;
    std::cout<<"  phiID: "<<phiID<<", thetaID "<<thetaID<<", rhoID "<<rhoID<<", cellID "<<cellID<<std::endl;
    std::cout<<"  Cell position: ("<<HitCellPos.x()<<", "<<HitCellPos.y()<<", "<<HitCellPos.z()<<std::endl;
    std::cout<<" (theta, phi, rho) = "<<"("<<HitCellPos.theta()<<", "<<HitCellPos.phi()<<", "<<HitCellPos.mag()<<") "<<std::endl;
    std::cout<<"  Neighbor cell count: "<<cellIDvec.size()<<std::endl;
    double responseSum = 0;
    for(size_t i=0; i<cellIDvec.size(); ++i) {
      phiID = m_segmentation->decoder()->get(cellIDvec[i], "phi");
      thetaID = m_segmentation->decoder()->get(cellIDvec[i], "theta");
      rhoID = m_segmentation->decoder()->get(cellIDvec[i], "rho");
      std::cout<<"  Neighbor cell "<<i<<": phiID "<<phiID<<", thetaID "<<thetaID<<", rhoID "<<rhoID<<", cellID "<<cellIDvec[i]<<", response "<<responseVec[i]<<std::endl;
      responseSum += responseVec[i];
    }
    std::cout<<"  Response sum: "<<responseSum<<std::endl;
#endif

    // Create the hits and accumulate contributions from multiple steps
    //
    Geant4HitCollection* coll = collection(m_collectionID);
    for (std::size_t i = 0; i < cellIDvec.size(); ++i) {
      const CellID hitCellID = cellIDvec[i];
      G4double longitudinalAttenuation = 1.;
      if (m_userData.fiberAttenuationLength > 0. && m_userData.useLightResponseFunction) {
        longitudinalAttenuation = std::exp(-(m_userData.outerRadius-global.perp()) / m_userData.fiberAttenuationLength);
      }
      // std::cout<<"  Response for cell "<<cellIDvec[i]<<": "<<responseVec[i]<<" / "<<longitudinalAttenuation<<std::endl;

      const G4double step_E = responseVec[i] * longitudinalAttenuation * aStep->GetTotalEnergyDeposit();
      Geant4Calorimeter::Hit* hit = coll->findByKey<Geant4Calorimeter::Hit>(hitCellID); // the hit

      if (!hit) { // if the hit does not exist yet, create it
        hit = new Geant4Calorimeter::Hit();
        hit->cellID = hitCellID;
        hit->position = cellPosVec[i]; // this should be assigned only once
        hit->energyDeposit = step_E;
        coll->add(hitCellID, hit); // add the hit to the hit collection
      } else {                 // if the hit exists already, increment its fields
        hit->energyDeposit += step_E;
      }

      // Add calo hit contributions
      // Note: Only add contribution for the central cell, to reduce the file size. 
      if(i==0){
        Geant4Calorimeter::Hit::Contribution contrib;
        contrib.trackID = aStep->GetTrack()->GetTrackID();
        contrib.pdgID = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
        contrib.deposit = aStep->GetTotalEnergyDeposit();
        // contrib.deposit = step_E; 
        contrib.time = aStep->GetPreStepPoint()->GetGlobalTime();
        contrib.x = global.x();
        contrib.y = global.y();
        contrib.z = global.z();
        hit->truth.emplace_back(contrib);
      }
    }

    return true;
  } // end of Geant4SensitiveAction::process() method specialization

} // namespace sim
} // namespace dd4hep

DECLARE_GEANT4SENSITIVE(GrainitaCaloSDAction)

//**************************************************************************
