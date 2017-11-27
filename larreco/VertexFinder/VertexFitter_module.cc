////////////////////////////////////////////////////////////////////////
// Class:       VertexFitter
// Plugin Type: producer (art v2_07_03)
// File:        VertexFitter_module.cc
//
// Author: Giuseppe Cerati, cerati@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "lardataobj/RecoBase/Vertex.h"

#include "larreco/RecoAlg/Geometric3DVertexFitter.h"

#include <memory>

namespace trkf {

  class VertexFitter : public art::EDProducer {
  public:

    struct Inputs {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> inputPFParticleLabel {
       Name("inputPFParticleLabel"),
       Comment("Label of recob::PFParticle Collection to be fit")
      };
      fhicl::Atom<art::InputTag> inputTracksLabel {
	Name("inputTracksLabel"),
	Comment("Label of recob::Track Collection associated to PFParticles")
      };
    };

    struct Config {
      using Name = fhicl::Name;
      fhicl::Table<VertexFitter::Inputs> inputs {
	Name("inputs"),
      };
      fhicl::Table<Geometric3DVertexFitter::Options> options {
	Name("options")
      };
      fhicl::Table<TrackStatePropagator::Config> propagator {
	Name("propagator")
      };
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit VertexFitter(Parameters const & p);
    
    // Plugins should not be copied or assigned.
    VertexFitter(VertexFitter const &) = delete;
    VertexFitter(VertexFitter &&) = delete;
    VertexFitter & operator = (VertexFitter const &) = delete;
    VertexFitter & operator = (VertexFitter &&) = delete;
    
    void produce(art::Event & e) override;

  private:
    art::InputTag pfParticleInputTag;
    art::InputTag trackInputTag;    
    Geometric3DVertexFitter fitter;
  };
}

trkf::VertexFitter::VertexFitter(Parameters const & p)
  : pfParticleInputTag(p().inputs().inputPFParticleLabel())
  , trackInputTag(p().inputs().inputTracksLabel())
  , fitter(p().options,p().propagator)
{
  produces<std::vector<recob::Vertex> >();
  produces<art::Assns<recob::PFParticle, recob::Vertex> >();
  produces<art::Assns<recob::Vertex, recob::Track, recob::TrackVertexMeta> >();
}

void trkf::VertexFitter::produce(art::Event & e)
{

  using namespace std;

  auto outputVertices = make_unique<vector<recob::Vertex> >();
  auto outputPFVxAssn = make_unique<art::Assns<recob::PFParticle, recob::Vertex> >();
  auto outputVxTkMtAssn = make_unique<art::Assns<recob::Vertex, recob::Track, recob::TrackVertexMeta> >();

  const auto& inputPFParticle = e.getValidHandle<vector<recob::PFParticle> >(pfParticleInputTag);
  auto assocTracks = unique_ptr<art::FindManyP<recob::Track> >(new art::FindManyP<recob::Track>(inputPFParticle, e, trackInputTag));

  // PtrMakers for Assns
  art::PtrMaker<recob::Vertex> vtxPtrMaker(e, *this);

  for (size_t iPF = 0; iPF < inputPFParticle->size(); ++iPF) {
    //
    art::Ptr<recob::PFParticle> pfp(inputPFParticle, iPF);
    if (pfp->IsPrimary()==false || pfp->NumDaughters()<2) continue;
    vector< art::Ptr<recob::Track> > tracks;
    auto& pfd = pfp->Daughters();
    for (auto ipfd : pfd) {
      vector< art::Ptr<recob::Track> > pftracks = assocTracks->at(ipfd);
      for (auto t : pftracks) {
	tracks.push_back(t);
      }
    }
    if (tracks.size()<2) continue;
    //
    VertexWrapper vtx = fitter.fitTracks(tracks);
    if (vtx.isValid()==false) continue;
    //
    auto meta = fitter.computeMeta(vtx, tracks);
    //
    // Fill the output collections
    //
    outputVertices->emplace_back(vtx.vertex());
    const art::Ptr<recob::Vertex> aptr = vtxPtrMaker(outputVertices->size()-1);
    outputPFVxAssn->addSingle( art::Ptr<recob::PFParticle>(inputPFParticle, iPF), aptr);
    //
    size_t itt = 0;
    for (auto t : tracks) {
      outputVxTkMtAssn->addSingle(aptr, t, meta[itt]);
      itt++;
    }
  }
  //
  e.put(std::move(outputVertices));
  e.put(std::move(outputPFVxAssn));
  e.put(std::move(outputVxTkMtAssn));

}

DEFINE_ART_MODULE(trkf::VertexFitter)
