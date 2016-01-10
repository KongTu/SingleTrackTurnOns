# SingleTrackTurnOns

This is a simple setup for checking the single track trigger efficiency in <strong> data </strong>. One needs to aware that
the trigger efficiency is combined HLT and L1 efficiency altogether. Also, the trigger efficiency is defined as 

<pre><code> leading pT distribution (higher threshold/lower threshold) </code></pre>

where the higher threshold cannot have any prescale, and the lowest threshold is calculated w.r.t MB trigger.

STEPS:
- Checkout at least 75X, and git clone this repo <pre><code> git clone https://github.com/KongTu/SingleTrackTurnOns.git </code></pre>
- Compile. <pre><code> scram b -j4 </code></pre>
- Go to /SingleTrackTurnOns/turnOnMaker/test, and replace the input file name, double check the track quality cuts and calo matching resolution parameter, before you do: <pre><code> cmsRun HLTandL1singletrack_crab_cfg.py </code></pre>
- One can use <strong> /macros/L1plusHLTSingleTrackTurnOns.C </strong> to plot efficiency as function of pT.



