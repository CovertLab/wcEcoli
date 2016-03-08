import os
from simulation import EcoliSimulation

class EcoliLiveSimulation(EcoliSimulation):
    def startup(self):
        # Perform initial mass calculations
        for state in self.states.itervalues():
            state.calculatePreEvolveStateMass()
            state.calculatePostEvolveStateMass()

        # Perform initial listener update
        for listener in self.listeners.itervalues():
            listener.initialUpdate()

        # Start logging
        for logger in self.loggers.itervalues():
            logger.initialize(self)
            
    def step(self, steps=1, time=None, sync_every=None):
        startTime = self._timeTotal
        stepCount = 0
        
        while True:
            stepCount += 1
            
            if steps > 0 and stepCount > steps:
                break
            
            if time is not None and startTime + time < self._timeTotal:
                break
            
            if self._cellCycleComplete:
                print "Cell cycle complete"
                break
                
            if sync_every is not None and stepCount % sync_every == 0:
                self._syncwrites()

            self._simulationStep += 1
            self._timeTotal += self._timeStepSec
            self._evolveState()
        
        self._syncwrites()
            
    def _syncwrites(self):
        # If python didn't want OO-breaking manipulation of private data,
        # it would have had visibility specification.

        for obj, saveFile in dict(self.loggers['Disk'].saveFiles, mainFile=self.loggers['Disk'].mainFile).iteritems():
            if saveFile._columns is not None:
                for name, column in saveFile._columns.iteritems():
                    column._data.flush()
                    column._offsets.flush()
                    os.fsync(column._data)
                    os.fsync(column._offsets)
    
    def shutdown(self):
        # Run post-simulation hooks
        for hook in self.hooks.itervalues():
            hook.finalize(self)

        # Divide mother into daughter cells
        self._divideCellFunction()

        # Finish logging
        for logger in self.loggers.itervalues():
            logger.finalize(self)