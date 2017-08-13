"""Helper object for files containing one or more MCParticle collections"""

from DDSim.Helper.Input import Input

class MCParticles( Input ):
  """Configuration for the generator-level InputFiles"""
  def __init__( self ):
    super(MCParticles, self).__init__()
    self._parameters["MCParticleCollectionName"] = "MCParticle"

  @property
  def mcParticleCollectionName( self ):
    """Set the name of the collection containing the MCParticle input.
    Default is "MCParticle".
    """
    return self._parameters["MCParticleCollectionName"]
  @particlesPerEvent.setter
  def particlesPerEvent( self, val ):
    self._parameters["MCParticleCollectionName"] = val

