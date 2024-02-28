//------------------------------------------
//	Analysis of long current injections
//	by Igor Delvendahl
//  and Philipp O'Neill
//  Subversion: For cultured human cells
//------------------------------------------

#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=CC_Analysis

static constant SPACER = 0.005


menu "Macros"
	"Analyse CC protocol", CC_Analysis#main()
end


static function main()
	// hard-coded parameters
	Variable TimeToEvaluate = 0.050 // time window used for input resistance (=50ms) and end of tau_membrane fit

	// default parameters
	String Wavestring
	Variable MINpeak = -0.01
	Variable dVdTthreshold = 2.5
	Variable dVdTthresholdFirstAP = 10
	Variable APsForFrequency = 5
	Variable SlidingAverageLevels = 1
	Variable SlidingAveragePeaks = 1
	Variable DetectionSpacer = 0.02
	String Protocol
	
	Prompt Wavestring, "Select first wave of series:", popup, WaveList("*", ";", "MINROWS:500") // show only data waves
	Prompt Protocol, "Protocol", popup, "Increase;Decrease;Toggle;"
	Prompt MINpeak, "AP peak minimum (V):"
	Prompt APsForFrequency, "Number of APs for frequency analysis:"
	Prompt dVdTthreshold, "AP threshold (V/s):"
	Prompt dVdTthresholdFirstAP, "AP threshold first AP (V/s):"
	// Set the box size for the sliding average for FindLevel and FindPeak functions.
	Prompt SlidingAverageLevels, "Sliding Average Thresh (levels):"
	Prompt SlidingAveragePeaks, "Sliding Average Thresh (peaks):"
	Prompt DetectionSpacer, "Spacer for AP detection"
	DoPrompt "CC analysis input parameters" Wavestring, Protocol, MINpeak, APsForFrequency, dVdTthreshold, dVdTthresholdFirstAP, SlidingAverageLevels, SlidingAveragePeaks, DetectionSpacer
	
	if(V_Flag)
		Print "User quit!"
		return -1								// User canceled	
	elseif( !IsHekaWave(Wavestring) )			// check if HEKA wave
		Print "Not a HEKA data wave."
		return -2
	endif
	
	Variable/G MINpeak_g = MINpeak
	Variable/G dVdTthreshold_g = dVdTthreshold
	Variable/G SlidingAverageLevels_g = SlidingAverageLevels
	Variable/G SlidingAveragePeaks_g = SlidingAveragePeaks

	
	String ListOfWaves = GetWavesFromHekaSeries(wavestring)  //make a list of all wave names from a given HEKA series
	variable NumberOfWaves = Itemsinlist(ListOfWaves)
	
	Make/O/T/N=(NumberOfWaves) Waves = StringFromList(p, ListOfWaves)

	String/G Filename = GetFilename(Waves[0])
	print Filename
	
	//get start and duration of current injection from first wave
	Variable StartOfCurrInject = GetSegmentDuration(Waves[0], 1)
	Variable DurOfCurrInject = GetSegmentDuration(Waves[0], 2)
	variable EndOfCurrInject = StartOfCurrInject + DurOfCurrInject

	// make current protocol wave
	Make/O/N=(NumberOfWaves) InjectedCurrent
	Variable currentstep
	if(NumberOfWaves == 1)	
		currentstep = GetCurrentIncreaseStep(Waves[0])
		InjectedCurrent = currentstep
	else 
		strswitch(Protocol)	
			case "Toggle":		// Toggle protocol (first sweep 0, then alternate positive & negative)
				currentstep = GetCurrentIncreaseStep(Waves[1])
				InjectedCurrent = MakeCurrentProtocolWave(NumberOfWaves, CurrentStep)
				break						
			default:		// increase or decrease: exact current injections are stored in wavenote & can be extracted directly
				currentstep = GetCurrentIncreaseStep(Waves[1]) - GetCurrentIncreaseStep(Waves[0])
				InjectedCurrent = GetCurrentIncreaseStep(Waves[p])
		endswitch
	endif
	
	// generate output waves
	make/D/O/N=(NumberOfWaves) V_rest=NaN, V_mean=NaN, deltaV=NaN, R_in=NaN,tau_membrane=NaN, SagRatio=NaN // features of negative current injections
	make/D/O/N=(NumberOfWaves) APThreshold = NaN, APamplitude = NaN, APduration = NaN // single AP features
	make/D/O/N=(NumberOfWaves) APnumber=NaN, maxAPfrequency=NaN, avgAPfrequency=NaN, adaptIndex=NaN, attenuation=NaN, broadening=NaN, isi_cv=NaN, pause=NaN, delay=NaN, APlatency=NaN, dVdt_ratio=NaN, lowpoint_lastAP = NaN, lowpoint_firstAP = NaN
	make/D/O/N=(NumberOfWaves) IFCs = NaN, ISImin = NaN, ISImax = NaN, ISImean = NaN, ISImedian = NaN, LastPeak = NaN
	make/D/O/N=0 spikesX, spikesY

	variable j
	display
	for(j = 0; j < NumberOfWaves; j += 1)
		wave voltageWave = $(Waves[j])
		printf "Analysed wave: %s\r" nameofwave(voltageWave)
		
		V_rest[j] = mean(voltageWave, SPACER, (StartOfCurrInject - SPACER)) * 1e3
		V_mean[j] = mean(voltageWave, (EndOfCurrInject - TimeToEvaluate - SPACER), (EndOfCurrInject - SPACER)) * 1e3
		deltaV[j] = V_mean[j] - V_rest[j]
		R_in[j] = deltaV[j] / InjectedCurrent[j] * 1e3

		// generate result wave for single AP features
		String ThresholdCrossingsref = "ThresholdCrossings" + num2str(j)
		String PeakValuesref = "PeakValues" + num2str(j)
		String PeakTimesref = "PeakTimes" + num2str(j)
		String APthresholdref = "APthreshold_all" + num2str(j)
		String APthresholdpositionref = "APthresholdposition" + num2str(j)
		String APamplituderef = "APamplitude_all" + num2str(j)
		String APdurationref = "APduration_all"  + num2str(j)
		String hfmaxleftref = "hfmaxleft" + num2str(j)
		String hfmaxrightref = "hfmaxright" + num2str(j)
		String lowpointref = "lowpoint_all" + num2str(j)
		String lowpointpositionref  = "lowpointposition" + num2str(j)
	
	
		// Smoothing for AP detection and differentiating for slope-wave dVdt
		Wave smoothedAPWave
		Duplicate/FREE/O voltageWave smoothedAPWave		
		Differentiate smoothedAPWave /D=dVdt
		Differentiate dVdt /D=secondDeriv
		
		
		if(InjectedCurrent[j] >= 0)		//positive current injection = potential APs
			
			// Treshold crossings; holds the X-values where MINpeak was crossed
			wave ThresholdCrossings
			FindLevels/B=(SlidingAverageLevels)/EDGE=1/DEST=ThresholdCrossings/Q/R=(StartOfCurrInject, EndOfCurrInject) voltageWave, MINpeak
			
			// If there are potential APs: find all relevant values.
			if (numpnts(ThresholdCrossings) > 0)
				
				// Find crossing of MINPeak
				make/D/O/N=(numpnts(ThresholdCrossings)) PeakTimes = NaN, PeakValues = NaN
				variable i
				for (i=0; i<numpnts(ThresholdCrossings);i+=1) // Find all PeakTimes in time coordinates as well as the peak values
					if (i > numpnts(ThresholdCrossings))
						FindPeak/Q/P/B=(SlidingAveragePeaks)/R=(ThresholdCrossings[i], ThresholdCrossings[i+1])/M=(MINpeak) voltageWave
					else
						FindPeak/Q/P/B=(SlidingAveragePeaks)/R=(ThresholdCrossings[i], EndOfCurrInject)/M=(MINpeak) voltageWave
					endif
				
					if (V_flag == 0)
						PeakValues[i] = V_PeakVal
						PeakTimes[i] = pnt2x(voltageWave, V_PeakLoc)
					else			
						PeakTimes[i] = NaN
						PeakValues[i] = NaN
						ThresholdCrossings[i] = NaN
					endif
				
				endfor
				
				// When no peak could be assigned to a threshold crossing: delete entry
				WaveTransform/O zapNaNs ThresholdCrossings
				WaveTransform/O zapNaNs PeakValues
				WaveTransform/O zapNaNs PeakTimes	
			
				// Remove all entries that do not have a slope steeper than 20V/s (or 10V/s, ...), and define thresholds if found
				make/D/O/N=(numpnts(ThresholdCrossings)) APthreshold_all=NaN, APthresholdposition = NaN
				for(i=0; i < numpnts(ThresholdCrossings); i+=1)		
					if (i != 0)
						FindLevel/Q/P/EDGE=1/B=(SlidingAverageLevels)/R=(PeakTimes[i], PeakTimes[i]-DetectionSpacer) dVdt, dVdTthreshold
					else
						FindLevel/Q/P/EDGE=1/B=(SlidingAverageLevels)/R=(PeakTimes[i], PeakTimes[i]-DetectionSpacer) dVdt, dVdTthresholdFirstAP
					endif
		
					if (V_flag == 0)
						APthreshold_all[i] = voltageWave[V_LevelX]
						APthresholdposition[i] = pnt2x(voltageWave, V_LevelX)
					else			
						PeakTimes[i] = NaN
						PeakValues[i] = NaN
						ThresholdCrossings[i] = NaN
						APthreshold_all[i] = NaN
						APthresholdposition[i] = NaN
					endif
				endfor
			
				WaveTransform/O zapNaNs ThresholdCrossings
				WaveTransform/O zapNaNs PeakValues
				WaveTransform/O zapNaNs PeakTimes
				WaveTransform/O zapNaNs APthreshold_all
				WaveTransform/O zapNaNs APthresholdposition
			
				make/D/O/N=(numpnts(ThresholdCrossings)) APamplitude_all=NaN, APduration_all=NaN, hfmaxleft = NaN, hfmaxright = NaN, lowpoint_all = NaN, lowpointposition = NaN
				for(i=0; i < numpnts(ThresholdCrossings); i+=1)
				
					// Find AP amplitude
					APamplitude_all[i] = PeakValues[i] - APthreshold_all[i]
					variable halfamplitude = APthreshold_all[i] + APamplitude_all[i] / 2
					
					// Find halfamplitude in the rising phase
					FindLevel/B=(SlidingAverageLevels)/EDGE=1/Q/R=(APthresholdposition[i], PeakTimes[i]) voltageWave, halfamplitude
					hfmaxleft[i] = V_LevelX
				
					// Find halfamplitude in the decay phase
					FindLevel/B=(SlidingAverageLevels)/EDGE=2/Q/R=(PeakTimes[i], PeakTimes[i]+DetectionSpacer) voltageWave, halfamplitude
					hfmaxright[i] = V_LevelX
				
					// Find FWHM
					APduration_all[i] =  (hfmaxright[i] - hfmaxleft[i]) * 1e6
					
					// Find lowest point in hyperpol phase
					WaveStats/Q/R=(PeakTimes[i], PeakTimes[i]+0.05) voltageWave
					lowpoint_all[i] = V_min
					lowpointposition[i] = V_minloc
				endfor
				
				// Display voltage trace with features used for analysis marked
				duplicate ThresholdCrossings, $ThresholdCrossingsref
				duplicate PeakValues, $PeakValuesref
				duplicate PeakTimes, $PeakTimesref
				duplicate APthreshold_all, $APthresholdref
				duplicate APthresholdposition, $APthresholdpositionref
				duplicate APamplitude_all, $APamplituderef
				duplicate APduration_all, $APdurationref
				duplicate hfmaxleft, $hfmaxleftref
				duplicate hfmaxright, $hfmaxrightref
				duplicate lowpoint_all, $lowpointref
				duplicate lowpointposition, $lowpointpositionref
				
				String ydat2ref = "ydat2_" + num2str(j)
				String ydat3ref = "ydat3_" + num2str(j)
				String ydat4ref = "ydat4_" + num2str(j)
				String ydat5ref = "ydat5_" + num2str(j)
				String ydat6ref = "ydat6_" + num2str(j)
				
				//display graph with voltage traces
				variable NumOfTraces = numpnts(Waves)
				AppendToGraph $(Waves[j])
				ModifyGraph rgb($(Waves[j]))=(1,34817,52428),btLen=4,gFont="Verdana",nticks(left)=3
				SetAxis/A/N=1/E=0 left
				
				if (numpnts(Peaktimes))
					// append the threshold crossings, Peak times etc. to graph to visualise detection
					make/O/D/N=(numpnts('PeakTimes')) ydat2
					ydat2 = VoltageWave('PeakTimes')
					duplicate ydat2 $ydat2ref
					AppendToGraph $ydat2ref vs $PeakTimesref
					ModifyGraph mode($ydat2ref)=3,marker($ydat2ref)=19, msize($ydat2ref)=1.5
					
					make/O/D/N=(numpnts('APthresholdposition')) ydat3
					ydat3 = voltageWave('APthresholdposition')
					duplicate ydat3 $ydat3ref
					AppendToGraph $ydat3ref vs $APthresholdpositionref
					ModifyGraph mode($ydat3ref)=3,marker($ydat3ref)=19, msize($ydat3ref)=1.5, rgb($ydat3ref)=(0,39168,0)
					
					make/O/D/N=(numpnts('hfmaxleft')) ydat4
					ydat4 = voltageWave('hfmaxleft')
					duplicate ydat4 $ydat4ref
					AppendToGraph $ydat4ref vs $hfmaxleftref
					ModifyGraph mode($ydat4ref)=3,marker($ydat4ref)=19, msize($ydat4ref)=1.3, rgb($ydat4ref)=(0,0,0)
					
					make/O/D/N=(numpnts('hfmaxright')) ydat5
					ydat5 = voltageWave('hfmaxright')
					duplicate ydat5 $ydat5ref
					AppendToGraph $ydat5ref vs $hfmaxrightref
					ModifyGraph mode($ydat5ref)=3,marker($ydat5ref)=19, msize($ydat5ref)=1.3, rgb($ydat5ref)=(65280,65280,0)
					
					make/O/D/N=(numpnts('lowpointposition')) ydat6
					ydat6 = voltageWave('lowpointposition')
					duplicate ydat6 $ydat6ref
					AppendToGraph $ydat6ref vs $lowpointpositionref
					ModifyGraph mode($ydat6ref)=3,marker($ydat6ref)=19, msize($ydat6ref)=1.3, rgb($ydat6ref)=(63232,43008,61184)
				endif
			endif // all APs and associated values found and appended to graph: Save values we will need for the results
			
			variable NumberOfEvents = numpnts(PeakTimes)
			
			if(NumberOfEvents > 0)
				// Get the last peaktime (for deciding which last hyperpol to take later)
				LastPeak[j] = PeakTimes[numpnts(PeakTimes)-1]
				
				// get number of APs and avg. frequency
				APnumber[j] = numpnts(ThresholdCrossings)
				avgAPfrequency[j] = APnumber[j] / DurOfCurrInject
				
				// get maxAPfrequency
				if (APnumber[j] == 0 || APnumber[j] == 1)
					maxAPfrequency[j] = avgAPfrequency[j]
				elseif (APnumber[j] == 2 || APnumber[j] == 3 ||APnumber[j] == 4)
					variable lastpoint = APnumber[j] - 1
					maxAPfrequency[j] = APnumber[j] / (PeakTimes[lastpoint] - PeakTimes[0])
				else
					maxAPfrequency[j] = APsForFrequency / (PeakTimes[APsForFrequency-1] - PeakTimes[0])
				endif
				
				// Get dVdt ratio, frequencies and latency
				dVdt_ratio[j] = abs(wavemax(dVdt, StartOfCurrInject, inf) / wavemin(dVdt, StartOfCurrInject, inf))
				APlatency[j] = PeakTimes[0] - StartOfCurrInject
				APthreshold[j] = APthreshold_all[0]													// read out the first AP threshold in each sweep
				APamplitude[j] = APamplitude_all[0]													// read out the first AP amplitude in each sweep
				APduration[j] = APduration_all[0]
				
				// Get last undershoot
				if (PeakTimes[numpnts(Peaktimes)-1]+0.004 <= EndOfCurrInject) //for last undershoot: sometimes will fall outside of the current injection... if the last AP too close to end of CI, take undershoot from the second last AP
					lowpoint_lastAP[j] = lowpoint_all[numpnts(lowpoint_all)-1]  * 1e3
				else
					lowpoint_lastAP[j] = lowpoint_all[numpnts(lowpoint_all)-2]  * 1e3
				endif
				
				// Get first undershoot
				lowpoint_firstAP[j] = lowpoint_all[0]  * 1e3
				
				// get adaptation index, attenuation, broadening, ISI CV, pause and delay...
				if(NumberOfEvents > 1)
					// Generate spikeY & spikeX for plot later
					insertpoints/M=0 numpnts(spikesY), APnumber[j], spikesY
					spikesY[numpnts(spikesY) - APnumber[j], numpnts(spikesY) - 1] = InjectedCurrent[j]
					concatenate/NP {ThresholdCrossings}, spikesX
					
					wave ISIs =  GetSpikeISIs(PeakTimes)
					WaveStats/Q ISIs
					ISImin[j] = V_min
					ISImax[j] = V_max
					ISImean[j] = V_avg
					ISImedian[j] = StatsMedian(ISIs)
					adaptIndex[j] = GetAdaptationIndex(ISIs)											// Adaptation index
					attenuation[j] = 1-(APamplitude_all[NumberOfEvents-1]/APamplitude_all[0])	// Attenuation
					broadening[j] = APduration_all[NumberOfEvents-1]/APduration_all[0]			// Broadening
					isi_cv[j] = GetISICV(ISIs)															// ISI Coefficient of variation
					pause[j] = GetSpikePause(ISIs)														// spike pause?
					delay[j] = (APlatency[j] > mean(ISIs)) ? 1 : 0									// delayed onset of spikes?
					IFCs[j] = GetIFC(PeakTimes, EndOfCurrInject, StartOfCurrInject)			// comparison last to first 500us
				else
					insertpoints/M=0 numpnts(spikesX), 1, spikesX, spikesY
					spikesY[numpnts(spikesY) - 1] = InjectedCurrent[j]
					spikesX[numpnts(spikesY) - 1] = NaN
				endif
			
			else // Sweep without any APs:
				LastPeak[j] = 0
				APnumber[j] = 0
				avgAPfrequency[j] = 0
				maxAPfrequency[j] = 0				
				dVdt_ratio[j] = 0
				APlatency[j] = 0
				APthreshold[j] = 0												
				APamplitude[j] = 0										
				APduration[j] = 0
				lowpoint_lastAP[j] = 0
				lowpoint_firstAP[j] = 0
				ISImin[j] = 0
				ISImax[j] = 0
				ISImean[j] = 0
				ISImedian[j] = 0
				adaptIndex[j] = 0											
				attenuation[j] = 0	
				broadening[j] = 0			
				isi_cv[j] = 0													
				pause[j] = 0														
				delay[j] = 0								
				IFCs[j] = 0			

			endif
		
		elseif(InjectedCurrent[j]  < 0)		//negative current injection = fit tau membrane and get sag ratio
			SagRatio[j] = GetSagRatio(voltageWave)
			APnumber[j] = 0
			avgAPfrequency[j] = 0
			try
				CurveFit/NTHR=0/Q/W=2/N exp_XOffset voltageWave [x2pnt(voltageWave, StartOfCurrInject), x2pnt(voltageWave, EndOfCurrInject - 0.05)] /D; AbortonRTE
				tau_membrane[j] = K2 * 1e3
			catch
				if (V_AbortCode == -4)
					Variable CurveFiterror = GetRTError(1)	// 1 to clear the error
					Printf "Error during curve fit: &s", GetErrMessage(CurveFiterror)
				endif
			endtry
			AppendToGraph $(Waves[j])
			ModifyGraph rgb($(Waves[j]))=(1,34817,52428),btLen=4,gFont="Verdana",nticks(left)=3
			SetAxis/A/N=1/E=0 left
		endif		
	endfor
	
	Sort InjectedCurrent,Waves,InjectedCurrent, V_rest, V_mean, deltaV, R_in,tau_membrane, SagRatio, APnumber, maxAPfrequency, avgAPfrequency, adaptIndex, attenuation, broadening, isi_cv, pause, delay, APlatency, APthreshold, APamplitude, APduration, dVdt_ratio, IFCs, ISImin, ISImax, ISImean, ISImedian, LastPeak, lowpoint_lastAP, lowpoint_firstAP

	// plotting:
	// display current-AP delay graph		
	Plot(APlatency, InjectedCurrent, ylab="AP latency (s)", fromZero=1)
		
	//display current-maxfrequency graph #1
	Plot(maxAPfrequency, InjectedCurrent, ylab="Maximum AP frequency (Hz)", fromZero=1)
		
	//display current-avgfrequency graph #2
	Plot(avgAPfrequency, InjectedCurrent, ylab="Average AP frequency (Hz)", fromZero=0)
	
	
	// Linear fit for f_i_slope	
	findlevel/EDGE=1/P/Q APnumber, 1
	if(V_flag)
		variable f_i_slope = 0
	else
		variable rheobase = InjectedCurrent[ceil(V_LevelX)] // find rheobase	

		variable LineFit_startP = BinarySearchInterp(InjectedCurrent, rheobase)-1
	
		Findvalue/T=1e-9/V=(wavemax(APnumber)) APnumber
		variable LineFit_endP = V_value
	
		CurveFit/NTHR=0/Q/W=2/N line avgAPfrequency[LineFit_startP, LineFit_endP] /X=InjectedCurrent /D 
		f_i_slope = K1
	endif
	
	
	//display current-AP number graph
	Plot(APnumber, InjectedCurrent, ylab="Number of APs", fromZero=1)
		
	//display R_in graph
	Plot(R_in, InjectedCurrent, ylab="Input resistance (MOhm)")		
	
	//display current-tau membrane graph
	Plot(tau_membrane, InjectedCurrent, ylab="Membrane time constant (ms)", fromZero=-1) 
		
	//display AP duration vs current
	Plot(APduration, InjectedCurrent, ylab="First AP FWHM", fromZero=1)
			
	//display current-voltage graph
	Plot(deltaV, InjectedCurrent, ylab="Delta voltage (mV)")	
	SetAxis/A/E=0 left
	//interpolate deltaV vs current
	Interpolate2/T=3/N=1000/F=0.01/Y=deltaV_interpolated InjectedCurrent, deltaV
	Differentiate deltaV_interpolated/D=deltaV_interpolated_DIF
	AppendToGraph deltaV_interpolated
	ModifyGraph rgb(deltaV_interpolated)=(0,0,0)
	
	//Display attenuation graph
	Plot(Attenuation, InjectedCurrent, ylab="Attenuation", fromZero=1)
	
	//Display ISI CV graph
	Plot(isi_CV, InjectedCurrent, ylab="ISI CV", fromZero=1)
	
	//display raster plot
	Display spikesY vs spikesX	
	ModifyGraph mode=3,marker=10,mrkThick=1,rgb=(1,34817,52428),btLen=4,gFont="Verdana",nticks(left)=3
	Label left "injected current (pA)"; Label bottom "s"; SetAxis bottom StartOfCurrInject,EndOfCurrInject; SetAxis/A/N=1/E=1 left

	//sort graph windows
	execute/Q/Z "Tilewindows/A=(2,6)/O=1/G=4/W=(10,30,1140,380)"

	//display table
	edit/W=(10,410,1140,810)/N=CompleteTable Waves,InjectedCurrent, V_rest, V_mean, deltaV, R_in,tau_membrane, SagRatio, APnumber, maxAPfrequency, avgAPfrequency, adaptIndex, attenuation, broadening, isi_cv, pause, delay, APlatency, APthreshold, APamplitude, APduration, dVdt_ratio, lowpoint_lastAP, lowpoint_firstAP, IFCs, ISImin, ISImax, ISImean, ISImedian
		
	//print results
	variable zero = BinarySearch(InjectedCurrent, 0)
	printf "\r---------------------------------\r"
	Printf "Processed %g waves \r", NumberOfWaves
	Printf "Start of current injection: %g ms \r", StartOfCurrInject * 1e3
	Printf "Duration of current injection: %g ms \r",(EndOfCurrInject - StartOfCurrInject) * 1e3
	Printf "Current increase step: %g pA \r", CurrentStep
	Printf "Input resistance from -%2d pA step: %g MOhm \r",CurrentStep, R_in[zero - 1] 
	Printf "Membrane time constant from -%2d pA step: %g ms \r", CurrentStep, tau_membrane[zero - 1]
	Printf "Rectification index: %g \r", R_in[zero + 1]  / R_in[zero-1] 
	Printf "Rheobase: %g pA \r", rheobase
	Printf "F-I Slope: %g Hz/pA \r", f_i_slope
	Printf "Maximum AP frequency: %g Hz \r", wavemax(maxAPfrequency)
	printf "\r---------------------------------\r"	
	
	//create output table and save/delete window
	saveresults(rheobase=rheobase, slope=f_i_slope)

end



// extracts source filename from wavenote
static function/S GetFilename(TheWaveString)
	string TheWaveString
	
	string some_wave = StringFromList(0, TheWaveString)
	string FilenameStr = greplist(note($some_wave), "(?s)(.*?)Filename|PMDatFile(.*?)")
	splitstring /E="([[:digit:]]+)-([[:digit:]]+)-([[:digit:]]+)_([[:digit:]]+)" FilenameStr

	return S_value
end


// returns the duration of the given segment
static function GetSegmentDuration(wavestring, segmentNo)
	string wavestring
	variable segmentNo
	
	string some_wave = StringFromList(0, wavestring)
	variable scalefactor = RunningWindows() ? 1 : 1e-3
	string SegmentName = "D" + num2str(segmentNo)
	string separatorStr = SelectString( RunningWindows(), "\r", ";" )
	variable SegmentDuration = str2num(stringbykey(SegmentName, note($some_wave), ":", separatorStr))

	return SegmentDuration * scalefactor
end


// returns the current increase step of segment 1
static function GetCurrentIncreaseStep(wavestring)
	string wavestring
	
	string some_wave = StringFromList(0, Wavestring)
	if(!waveexists($some_wave))
		print "wave", some_wave, "does not exist"
		return NaN
	endif
	variable scalefactor = RunningWindows() ? 1e3 : 1
	string WaveMetaData = note($some_wave)
	string separatorStr = SelectString( RunningWindows(), "\r", ";" )
	string V1 = stringbykey("V1", WaveMetaData, ":", separatorStr)
	string V2 = stringbykey("V2", WaveMetaData, ":", separatorStr)
	variable current_increase = str2num(V2) - str2num(V1)

	return current_increase * scalefactor
end


// generate wave with current injections (assumes toggle protocol)
static function MakeCurrentProtocolWave(NumberOfWaves, CurrentStep)
	variable NumberOfWaves
	variable CurrentStep
	
	variable i, j
	Make/O/N=(NumberOfWaves)/FREE TempCurrent
	TempCurrent[0] = 0
	for (i = 1; i < NumberOfWaves; i += 1)
		if (mod(i, 2) == 0)
			TempCurrent[i] = j * CurrentStep * (-1)
		else
			j += 1
			TempCurrent[i] = j * CurrentStep
		endif
	endfor
	
	return TempCurrent
end


// returns wave with ISIs (calculated from SpikeTimes)
Static Function/WAVE GetSpikeISIs(SpikeTimes)
	wave SpikeTimes
	
	variable NumOfSpikes = numpnts(SpikeTimes)
	Make/O/D/FREE/N=(NumOfSpikes - 1) ISI 	
	ISI = SpikeTimes[p + 1] - SpikeTimes[p]
	
	return ISI
End


// returns adaptation index
Static Function GetAdaptationIndex(ISI)
	wave ISI

	variable adaptation = 0
	variable i
	variable lenISI = numpnts(ISI)
	for(i = 0; i < (lenISI - 1); i += 1)
		adaptation += (ISI[i + 1] - ISI[i]) / (ISI[i + 1] + ISI[i])
	endfor	
	adaptation /=  (lenISI - 1)
	
	return adaptation
End


// returns CV of spike ISIs
Static Function GetISICV(ISI)
	wave ISI

	wavestats/Q ISI
	variable CV = v_sdev / v_avg
	
	return CV
End	


// check if spikes have a pause (= ISI > 3*ISI_pre and ISI_post)
Static Function GetSpikePause(ISI)
	wave ISI
	
	variable numISI = numpnts(ISI)
	if(numISI > 3)
		variable i = 1
		do
			if( ISI[i] > (3 * ISI[i - 1]) && ISI[i] > (3 * ISI[i + 1]) )
				return 1
			endif
			i += 1	
		while( i < (numISI -1) )
	endif
	
	return 0
End


// calculate intrinsic frequency change as defined by Masoli et al. 2020
Static Function GetIFC(PeakTime, EndTime, StartTime)
	wave PeakTime
	variable EndTime, StartTime
	
	variable TotalNumber = numpnts(PeakTime)
	variable FreqBegin = 0 
	variable FreqEnd = 0
	variable IFC
	variable i
	for (i = 0; i < TotalNumber; i += 1)
		if (PeakTime[i] >= StartTime && PeakTime[i] <= (StartTime + 0.5))
			FreqBegin += 1
		elseif (PeakTime[i] <= EndTime && PeakTime[i] >= (EndTime - 0.5))
			FreqEnd += 1
		endif
	endfor
	IFC = FreqEnd/FreqBegin
	
	return IFC
End


// returns sag ratio
Static Function GetSagRatio(InputWave)
	wave InputWave
	
	string ThisWave = nameofwave(InputWave)
	
	variable SagRatio
	Variable StartOfCurrInject = GetSegmentDuration(ThisWave,1)
	Variable DurOfCurrInject = GetSegmentDuration(ThisWave,2)

	Duplicate/O/FREE InputWave SmoothedWave
	Smooth 50, SmoothedWave

	Variable baseline = mean(SmoothedWave,0,StartOfCurrInject)
	variable HalfOfCurrentInjection = (StartOfCurrInject + DurOfCurrInject) / 2
	wavestats/M=1/Q/R=(StartOfCurrInject,StartOfCurrInject + HalfOfCurrentInjection) SmoothedWave
	Variable peak = V_min
	Variable steadystate = mean(SmoothedWave, (StartOfCurrInject + DurOfCurrInject) - 0.05, StartOfCurrInject + DurOfCurrInject)
	
	SagRatio = ((baseline-peak) - (baseline-steadystate)) / (baseline-peak)
	
	return SagRatio
End


// returns 0 if Macintosh, 1 if Windows
static function RunningWindows()
	string platform = UpperStr(igorinfo(2))
	
	return strsearch(platform, "WINDOWS", 0) >= 0
end


//make a list of all wave names from given HEKA series
static function/S GetWavesFromHekaSeries(wavestring)	
	string wavestring
	
	string prefix, mode
	variable pm_expt, pm_series, pm_sweep, pm_trace
	sscanf wavestring, "%[^_]%*[_]%d%*[_]%d%*[_]%d%*[_]%d%*[_]%s", prefix, pm_expt, pm_series, pm_sweep, pm_trace, mode
	string list = ""
	string wavetoadd
	string formatStr = SelectString( RunningWindows(), "%s_%1.0f_%1.0f_%1.0f_%1.0f_%s", "%s_%1.0f_%1.0f_%03.0f_%1.0f_%s" )
	
	do
		sprintf wavetoadd, formatStr, prefix, pm_expt, pm_series, pm_sweep, pm_trace, mode		
		if(waveexists($wavetoadd) == 0)
			break
		endif
		list = addlistitem(wavetoadd, list, ";", 999)
		pm_sweep += 1
	while(1)

	list = ListMatch(list, "*" + mode + "*")
	
	return  list
end


// plotting function for scatterplot
Static Function Plot(yWave, xWave, [ylab, xlab, fromZero])
	wave yWave, xWave
	string ylab, xlab
	variable fromZero
	
	if( ParamIsDefault(xlab) )
		xlab = "injected current (pA)"
	endif
	
	display yWave vs xWave
	ModifyGraph mode=3, marker=8, rgb=(0,0,0), msize=4, rgb=(1,34817,52428), btLen=4, gFont="Verdana", nticks(left)=3
	Label left ylab
	Label bottom xlab
	SetAxis/A/N=1/E=1 left
	if(fromZero == 1)
		SetAxis bottom 0,*
	elseif(fromZero == -1)
		SetAxis bottom *,0	
	endif	
End


//save results to file
	static function saveresults([rheobase, slope])
	variable rheobase, slope
	wave InjectedCurrent, V_rest, R_in, R_in_spline, tau_membrane, SagRatio, APnumber,maxAPfrequency,AvgAPfrequency,adaptIndex,attenuation,broadening,isi_cv,pause,delay,APlatency,APthreshold,APamplitude,APduration,dVdt_ratio,IFCs, ISImin, ISImax, ISImean, ISImedian, lowpoint_firstAP, lowpoint_lastAP, PeakTimes
	variable numResults = 30
	Make/O/N=(numResults)/T labels
	Make/O/N=(numResults)/D results
	
	
	variable zero = BinarySearch(InjectedCurrent, 0)
	variable minus20 = BinarySearch(InjectedCurrent, -20)
	Findvalue/T=1e-9/V=(wavemax(APnumber)) APnumber
	variable maxAPtrace = V_value
	
	
	SVAR Filename = :Filename
	variable cell_id
	string Filename2 = replacestring("-", Filename, "")
	Filename2 = replacestring("_", Filename2, "")
	sscanf Filename2, "%d", cell_id
	
	labels[0] = "cell_ID"; 							results[0] = cell_id
	labels[1] = "current step"; 						results[1] = InjectedCurrent[zero+1]
	labels[2] = "AP number"; 							results[2] = wavemax(APnumber)
	labels[3] = "max. AP frequency";					results[3] = wavemax(maxAPfrequency)
	labels[4] = "avg. AP frequency"; 				results[4] = wavemax(AvgAPfrequency)
	labels[5] = "adaptation";							results[5] = adaptIndex[maxAPtrace]
	labels[6] = "attenuation";						results[6] = attenuation[maxAPtrace]
	labels[7] = "broadening";							results[7] = broadening[maxAPtrace]
	labels[8] = "ISI CV";								results[8] = isi_cv[maxAPtrace]
	labels[9] = "Pause";								results[9] = pause[maxAPtrace]
	labels[10] = "Delay";								results[10] = delay[maxAPtrace]
	labels[11] = "AP latency";						results[11] = APlatency[maxAPtrace]
	labels[12] = "AP threshold";						results[12] = APthreshold[maxAPtrace] * 1e3
	labels[13] = "AP amplitude";						results[13] = APamplitude[maxAPtrace] * 1e3
	labels[14] = "AP fwhm";							results[14] = APduration[maxAPtrace]
	labels[15] = "AP dVdT ratio";						results[15] = dVdt_ratio[maxAPtrace]
	labels[16] = "IFC";									results[16] = IFCs[maxAPtrace]
	labels[17] = "Min ISI";							results[17] = ISImin[maxAPtrace]
	labels[18] = "Max ISI";							results[18] = ISImax[maxAPtrace]
	labels[19] = "Mean ISI";							results[19] = ISImean[maxAPtrace]
	labels[20] = "Median ISI";						results[20] = ISImedian[maxAPtrace]
	labels[21] = "First Undershoot";					results[21] = lowpoint_firstAP[maxAPtrace]
	labels[22] = "resting membrane potential"; 	results[22] = V_rest[zero]
	labels[23] = "input resistance"; 				results[23] = R_in[zero - 1]
	labels[24] = "Last Undershoot";					results[24] = lowpoint_lastAP[maxAPtrace]
	labels[25] = "membrane time constant";			results[25] = tau_membrane[zero - 1]
	labels[26] = "sag ratio";							results[26] = SagRatio[minus20]
	labels[27] = "rectification";						results[27] = R_in[zero + 1]  / R_in[zero-1] 
	labels[28] = "rheobase";							results[28] = rheobase
	labels[29] = "f-I slope";							results[29] = slope


	edit/N=ResultsTable labels, results
	
	Duplicate/O/FREE APnumber, APs_export
	Wavetransform zapNans APs_export
	
	string txt_results1, txt_results2
	wfprintf txt_results1, "%.8G\t"/R=[1,19], results
	wfprintf txt_results2, "%.8G\t", APs_export
	
	putscraptext Filename + "\t" + txt_results1 + "\t" + txt_results2
	PanelMake()
End


Static Function PanelMake()
	DoWindow CC_Protocol_SaveFiles
	if(V_Flag)
		Dowindow/F CC_Protocol_SaveFiles
		return 0
	endif
	
	NewPanel/W=(800,100,1000,290)/N=CC_Protocol_SaveFiles
	ModifyPanel/W=CC_Protocol_SaveFiles, frameStyle=3
	Button SaveResults title="Save Results Tables",pos={20,20},size={145,20},proc=CC_Analysis#SaveResultsTable, appearance = {os9}, font="Helvetica"
	Button SaveComplete title="Save Complete Table",pos={20,60},size={145,20},proc=CC_Analysis#SaveCompleteTable, appearance = {os9}, font="Helvetica"
	Button GetRseries title="Save Series Resistance",pos={20,100},size={145,20},proc=CC_Analysis#SaveRseries, appearance = {os9}, font="Helvetica"	
	Button cleanup title="Clean up",pos={20,140},size={145,20},proc=CC_Analysis#cleanup, appearance = {os9}, font="Helvetica"
End


// Save the Results Table
static function SaveResultsTable(buttonname)
	string buttonname
	SVAR Filename = :Filename
	string ActualName = Filename + "_ci_long_ResultsTable.csv"
	SaveTableCopy/I/W=ResultsTable/T=2 as ActualName
end


static function SaveCompleteTable(buttonname)
	string buttonname
	SVAR Filename = :Filename
	string ActualName = Filename + "_CCiv_CompleteTable.csv"
	SaveTableCopy/I/W=CompleteTable/T=2 as ActualName
end


static function SaveRseries(buttonname)
	string buttonname
	SVAR Filename = :Filename
	string name_of_Waves = WaveList("*", ";", "")
	string some_wave = StringFromList(0, name_of_Waves)
	variable Rseries = str2num(stringbykey("PMRs", note($some_wave)))*1e-6
	string ActualName = Filename + "_CCiv_Rseries.txt"
	make/D/O/N=1 Rser = Rseries
	save/O/J/I Rser as ActualName
end


// Variant of SaveGraphAsPicture() with name already prepared
Menu "Graph"
	"Save ci_long graph as...", ci_long_Analysis#SaveTheci_longGraph()
End


static function SaveTheci_longGraph()	
	SVAR Filename = :Filename
	string ActualName = Filename + "_Ci_Graph"
	variable fileformat = 4
	variable resolution = 600
	variable width = 10
	variable height = 4
	
	Prompt fileformat, "File format:",popup,("PDF;TIFF;JPEG;PNG;BMP;EPS")
	Prompt resolution, "Resolution (dpi):"
	Prompt width, "Image width (in):"
	Prompt height, "Image height (in):"
	DoPrompt "Save graph" fileformat,resolution,width,height
	if (V_Flag)
		Print "User quit!"
		return -1								// User canceled
	endif
	fileformat -= 9
	SavePICT/E=(fileformat)/RES=(resolution)/I/W=(0,0,width,height)/O as ActualName
End


// Kill the root data folder
Static function cleanup(buttonname)
	string buttonname
	KillGraphsandTables()
	KillDataFolder root:
End	


// close all Graph and Table windows
Static function KillGraphsAndTables()
	string WindowList = WinList("*", ";", "WIN:3")
	variable i
	variable numWindows = itemsinlist(WindowList)
	for(i=0; i<numWindows; i+=1)
		DoWindow/K $stringfromlist(i, WindowList)
	endfor
end


// check wave note if data was loaded from HEKA .dat file
Static Function IsHekaWave(WaveNamestring)
	string Wavenamestring

	if(grepstring(note($Wavenamestring), "v2x"))
		return 1
	else
		return 0
	endif
end