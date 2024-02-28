//------------------------------------------
//	Analysis of voltage step protocol
//	by Igor Delvendahl
//  
//------------------------------------------

#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName=IV

Menu "Macros"
	Submenu "Analyse IV"
		"Analyse vc IV",IV#Analyse()
		"Analyse recovery from inact",IV#AnalyseRecovery()
		"Analyse onset of inact",IV#AnalyseOnset()
		"Subtract traces",IV#SubtractTraces()
		"Analyse SS inactivation",IV#AnalyseInactivation()
		"analyse tail currents",IV#AnalyseTailCurrents()
	End
	Submenu "Analyse Ca current"
		"Analyse IV",IV#AnalyseCaIV()
		"Analyse inactivation",IV#Analyse_CaInactivation()
	End
End

//***********************************************************

static constant k_spacer = 0.0005
static strconstant ks_prefix = "PMPulse_1_"

//***********************************************************
static function Analyse()
	String ListOfWaves = SelectTraces(SeriesNum = GetSeriesNum())
	if(cmpstr(ListOfWaves,""))
		variable Segment = GetRelevantSegment(1)
		variable VoltageSegment =  GetVoltageSegment(1)
		MakeRsCorrection(ListOfWaves)
		EvaluateCurrents(ListOfWaves,Segment,VoltageSegment)
	endif
end
//***********************************************************
static function AnalyseRecovery()
	String ListOfWaves = SelectTraces(SeriesNum = GetSeriesNum())
	if(cmpstr(ListOfWaves,""))
		EvalRecovFromInact(ListOfWaves)
	endif
end
//***********************************************************
static function AnalyseOnset()
	String ListOfWaves = SelectTraces(SeriesNum = GetSeriesNum())
	if(cmpstr(ListOfWaves,""))
		EvalOnsetOfInact(ListOfWaves)
	endif
end
//***********************************************************
static function AnalyseInactivation()
	String ListOfWaves = SelectTraces(SeriesNum = GetSeriesNum())
	if(cmpstr(ListOfWaves,""))
		variable Segment = GetRelevantSegment(2)
		variable VoltageSegment =  GetVoltageSegment(1)
		EvaluateCurrents(ListOfWaves,Segment,VoltageSegment)
	endif
end

//***********************************************************
static function AnalyseTailCurrents()
	String ListOfWaves =SelectTraces(SeriesNum = GetSeriesNum())
	if(cmpstr(ListOfWaves,""))
		variable Segment = GetRelevantSegment(1)
		EvaluateTail(ListOfWaves,Segment)
	endif
end

//***********************************************************
static function AnalyseCaIV()
	String ListOfWaves = SelectTraces(SeriesNum = GetSeriesNum())
	if(cmpstr(ListOfWaves,""))
		variable Segment = GetRelevantSegment(1)
		MakeRsCorrection(ListOfWaves)
		EvaluateCalciumCurrents(ListOfWaves,Segment)
	endif
end

//***********************************************************
static function Analyse_CaInactivation()
	String ListOfWaves = SelectTraces(SeriesNum = GetSeriesNum())
	if(cmpstr(ListOfWaves,""))
		variable Segment = GetRelevantSegment(2)
		EvaluateDeactKinetics(ListOfWaves,Segment)
	endif
end

//***********************************************************
static function /S SelectTraces([SeriesNum])
	variable SeriesNum
	string seriesStr 
	
	if( ParamIsDefault(SeriesNum)  )
		seriesStr = "1"
	else
		seriesStr = num2str(SeriesNum)
	endif
	
	Prompt seriesStr, "series ("+ks_prefix+"X_1_1):"

	DoPrompt "IV analysis Input Parameters" seriesStr
	if (V_Flag)
		return ""
	endif
	KillLeakWaves()

	string ListOfWaves = WaveList(ks_prefix+seriesStr+"*",";","")
	if (!itemsinlist(ListOfWaves))
		DoAlert 0, "No traces found."
		return ""
	endif

	return ListOfWaves
end

//***********************************************************
static function GetRelevantSegment(defaultSegment)
	variable defaultSegment
	variable segment = defaultSegment
	Prompt segment, "PGF segment for analysis:"

	DoPrompt "Select Segment" segment
	if (V_Flag)
		return -1
	endif
	return segment
end
//***********************************************************
static function GetVoltageSegment(defaultSegment)
	variable defaultSegment
	variable segment = defaultSegment
	Prompt segment, "PGF segment for voltage:"

	DoPrompt "Select Segment" segment
	if (V_Flag)
		return -1
	endif
	return segment
end
//***********************************************************
static function EvaluateCurrents(ListOfWaves,RelevantSegment,VoltageSegment)
	string ListOfWaves
	variable RelevantSegment,VoltageSegment
	
	if(RelevantSegment < 0)
		DoAlert 0, "Non-existing segment"
		return -1
	endif
	string I_strings = ListMatch(ListOfWaves,"*Imon*")
	variable NumberOfWaves = Itemsinlist(I_strings)
	
	print GetMetaData(stringfromlist(0,I_strings))
	
	if(NumberOfWaves > 0)
		Make/O/T/N=(NumberOfWaves) I_traces = StringFromList(p,I_strings)
	endif

	variable StartOfProtocol,EndOfProtocol
	variable i
	for (i = 0; i < RelevantSegment; i += 1)
		StartOfProtocol += GetSegmentDuration(I_traces[1],i + 1)
	endfor
	EndOfProtocol = StartOfProtocol + GetSegmentDuration(I_traces[1],RelevantSegment + 1)

	Make/O/D/N=(NumberOfWaves) voltage,baseline, Na_peak, K_peak, K_ss,Na_Peak_2,Na_loc
	string StoreValue,NaPeakStr,KPeakStr,KssStr
	NaPeakStr = ""
	KPeakStr = ""
	KssStr = ""
	
	variable base, NegativeLevel,PeakVal
	string VoltageSegmentString = "V" + num2istr(VoltageSegment + 1)
	
	for(i = 0; i < NumberOfWaves; i += 1)
		wave raw_trace = $(I_traces[i])
		Duplicate/FREE raw_trace smth_trace
		Smooth 100, smth_trace

		voltage[i] = str2num(stringbykey(VoltageSegmentString,note($(I_traces[i])),":","\r"))
		base = mean(raw_trace,k_spacer,(StartOfProtocol - k_spacer))
		baseline[i] = base * 1e12

		NegativeLevel = base - (2e-11)
		variable peakPos = FindNegPeakPosition(smth_trace,StartOfProtocol+0.00025,StartOfProtocol + 0.05,NegativeLevel)
		variable peakAvg = mean(raw_trace, peakPos - 0.0001, PeakPos + 0.0001)
		Na_peak[i] = (peakAvg - base) * 1e12
		sprintf StoreValue, "%.8g\t", Na_peak[i]
		NaPeakStr +=  StoreValue
		Na_loc[i] = peakPos

		PeakVal = EvalPeak($(I_traces[i]),StartOfProtocol)
		Na_Peak_2[i] = PeakVal ? (PeakVal-base) * 1e12 : 0

		K_Peak[i] = (WaveMax(smth_trace,StartOfProtocol + 0.001,EndOfProtocol - 0.001) - base) * 1e12
		sprintf StoreValue, "%.8g\t", K_peak[i]
		KPeakStr +=  StoreValue

		K_ss[i] = (mean(smth_trace,(EndOfProtocol - 0.0102),(EndOfProtocol - 0.0002)) - base) * 1e12
		sprintf StoreValue, "%.8g\t", K_ss[i]
		KssStr +=  StoreValue
	endfor
	
	if (WindowExists("results") == FALSE)
		Edit/N=results/W=(0,475,600,675) voltage,baseline, Na_peak, K_peak, K_ss,Na_Peak_2
	endif	
	if (CheckIfWaveIsDisplayed(Na_loc) == FALSE)
		Display/K=1 Na_loc vs voltage
	endif
	if (CheckIfWaveIsDisplayed(Na_peak) == FALSE)	
		Display/K=1 Na_peak vs voltage
	endif
	if (CheckIfWaveIsDisplayed(K_peak) == FALSE)	
		Display/K=1 K_ss,K_peak vs voltage
	endif
	if (WindowExists("traces") == FALSE)
		Display/N=traces/K=1
		for(i = 0; i < NumberOfWaves; i += 1)
			AppendToGraph $(I_traces[i])
		endfor
	endif	
	
	string cmd = "TileWindows/A=(2,2)/O=1/W=(0,0,600,470)"
	execute cmd
	
	string clipboard = ""
	clipboard += NaPeakStr + "\t"
	clipboard += KPeakStr + "\t"
	clipboard += KssStr
	putscraptext clipboard
end
//***********************************************************

static function FindNegPeakPosition(inwave,startX,endX,level)	//functon to find inward current peak
	wave inwave
	variable startX,endX
	variable level
	
	variable prot_end = 0.15
	variable peak_max
	variable peak_loc
	
	variable search_start = startX + (k_spacer/10)
	variable windowsize = 0.0001
	do
		variable search_end = search_Start + windowsize
		FindPeak/B=20/M=(level)/N/Q/R=(search_start, search_end) inwave
		
		if(abs(V_PeakVal) > abs(peak_max))
			peak_loc = V_PeakLoc
			peak_max = V_PeakVal
		endif
		
		search_start = search_end
	while(search_start < prot_end)
	
	return peak_loc
end

//***********************************************************
static function EvalRecovFromInact(ListOfWaves)
	string ListOfWaves

	string I_strings = ListMatch(ListOfWaves,"*Imon*")
	variable NumberOfWaves = Itemsinlist(I_strings)
	
	print GetMetaData(stringfromlist(0,I_strings))
	
	if(NumberOfWaves > 0)
		Make/O/T/N=(NumberOfWaves) I_traces = StringFromList(p,I_strings)
	endif
		
	Make/O/D/N=(NumberOfWaves) baseline, Na_peak, Na_Peak_2,deltaT,Recovery
	string StoreValue,NaPeakStr
	NaPeakStr = ""
	
	variable i
	for(i = 0; i < NumberOfWaves; i += 1)	//loop through waves
		wave raw_trace = $(I_traces[i])
		Duplicate/FREE raw_trace smth_trace
		Smooth 1000, smth_trace
		
		variable pulse1_start = GetStartOfSegment(raw_trace,3)
		variable pulse1_end = GetEndOfSegment(raw_trace,3)
		variable pulse2_start = GetStartOfSegment(raw_trace,5)
		variable pulse2_end = GetEndOfSegment(raw_trace,5)
		
		deltaT[i] = pulse2_start - pulse1_end
				
		variable base = mean(raw_trace,k_spacer,(pulse1_start - k_spacer))
		baseline[i] = base * 1e12
		
		variable NegativeLevel = base - (3.5e-11)
		//pulse 1
		variable peakPos = FindNegPeakPosition(smth_trace,pulse1_start,pulse1_end,NegativeLevel)
		variable peakAvg = mean(raw_trace, peakPos - 0.0001, PeakPos + 0.0001)
		Na_peak[i] = (peakAvg - base) * 1e12
		
		//pulse 2
		variable base2 = mean(raw_trace, (pulse2_start - 0.0001), pulse2_start)
		peakPos = FindNegPeakPosition(smth_trace,pulse2_start,pulse2_end,NegativeLevel)
		peakAvg = mean(raw_trace, peakPos - 0.0001, PeakPos + 0.0001)
		variable peakval = (peakAvg - base2) * 1e12
		Na_peak_2[i] = (deltaT[i] > 0 && peakval <= 0) ? peakval : 0
		
		sprintf StoreValue, "%.8g\t", Na_peak_2[i]
		NaPeakStr +=  StoreValue		
	endfor
	
	variable Pulse1_avg = mean(Na_peak)
	Recovery = Na_Peak_2 / Pulse1_avg

	if (WindowExists("results") == FALSE)
		Edit/N=results/W=(0,475,600,725) deltaT,baseline,Na_peak,Na_Peak_2
	endif	
	if (CheckIfWaveIsDisplayed(Na_peak) == FALSE)	
		Display/K=1 Na_peak vs deltaT
		ModifyGraph mode=3,marker=8,msize=4,opaque=1
		ModifyGraph log(bottom)=1
	endif
	if (CheckIfWaveIsDisplayed(Recovery) == FALSE)	
		Display/K=1 Recovery vs deltaT
		ModifyGraph mode=4,marker=8,msize=4,opaque=1
	endif
	if (CheckIfWaveIsDisplayed(Na_Peak_2) == FALSE)	
		Display/K=1 Na_Peak_2 vs deltaT
		ModifyGraph mode=4,marker=8,msize=4,opaque=1
	endif
	if (WindowExists("traces") == FALSE)
		Display/N=traces/K=1
		for(i = 0; i < NumberOfWaves; i += 1)
			AppendToGraph $(I_traces[i])
		endfor
	endif	
	
	string cmd = "TileWindows/A=(2,2)/O=1/W=(0,0,600,450)"
	execute cmd
	
	string clipboard = ""
	NaPeakStr = RemoveListItem(0,NaPeakStr,"\t")	//remove amplitude of first sweep
	sprintf clipboard, "%.8g\t",Pulse1_avg			//avg amplitude of pulse 1
	clipboard += NaPeakStr							//combine for clipboard
	putscraptext clipboard
end

//***********************************************************
static function EvalOnsetOfInact(ListOfWaves)
	string ListOfWaves

	string I_strings = ListMatch(ListOfWaves,"*Imon*")
	variable NumberOfWaves = Itemsinlist(I_strings)
	
	print GetMetaData(stringfromlist(0,I_strings))
	
	if(NumberOfWaves > 0)
		Make/O/T/N=(NumberOfWaves) I_traces = StringFromList(p,I_strings)
	endif
		
	Make/O/D/N=(NumberOfWaves) baseline, Na_peak, deltaT, Onset
	string StoreValue,NaPeakStr
	NaPeakStr = ""
	
	variable i
	for(i = 0; i < NumberOfWaves; i += 1)	//loop through waves
		wave raw_trace = $(I_traces[i])
		Duplicate/FREE raw_trace smth_trace
		Smooth 1000, smth_trace

		variable pulse_start = GetStartOfSegment(raw_trace,4)
		variable pulse_end = GetEndOfSegment(raw_trace,4)
		
		deltaT[i] = pulse_start - GetEndOfSegment(raw_trace,2)
						
		variable base = mean(raw_trace, (GetEndOfSegment(raw_trace,2) - 0.0001), GetEndOfSegment(raw_trace,2))
		baseline[i] = base * 1e12
		variable NegativeLevel = base - (3.5e-11)

		variable peakPos = FindNegPeakPosition(smth_trace,pulse_start,pulse_end,NegativeLevel)
		variable peakAvg = mean(raw_trace, peakPos - 0.0001, PeakPos + 0.0001)
		variable peakval = (peakAvg - base) * 1e12
		Na_peak[i] = (peakval <= 0) ? peakval : 0
		
		sprintf StoreValue, "%.8g\t", Na_peak[i]
		NaPeakStr +=  StoreValue		
	endfor
	
	variable Pulse1_avg = mean(Na_peak)
	Onset = Na_Peak / Na_Peak[0]

	if (WindowExists("results") == FALSE)
		Edit/N=results/W=(0,475,600,725) deltaT,baseline,Na_peak
	endif	
	if (CheckIfWaveIsDisplayed(Na_peak) == FALSE)	
		Display/K=1 Na_peak vs deltaT
		ModifyGraph mode=3,marker=8,msize=4,opaque=1
		ModifyGraph log(bottom)=1
		
		Display/K=1 Na_Peak vs deltaT
		ModifyGraph mode=4,marker=8,msize=4,opaque=1
	endif
	if (CheckIfWaveIsDisplayed(Onset) == FALSE)	
		Display/K=1 Onset vs deltaT
		ModifyGraph mode=4,marker=8,msize=4,opaque=1
	endif
	if (WindowExists("traces") == FALSE)
		Display/N=traces/K=1
		for(i = 0; i < NumberOfWaves; i += 1)
			AppendToGraph $(I_traces[i])
		endfor
	endif	
	
	string cmd = "TileWindows/A=(2,2)/O=1/W=(0,0,600,450)"
	execute cmd

	putscraptext NaPeakStr
end

//***********************************************************
static function GetStartOfSegment(inwave,segment)
	wave inwave
	variable segment
		
	string wavestring = nameofwave(inwave)
	variable i,StartOfProtocol
	for (i = 1; i < segment; i += 1)
		StartOfProtocol += GetSegmentDuration(wavestring,i)
	endfor
	
	return StartOfProtocol
end
//***********************************************************
static function GetEndOfSegment(inwave,segment)
	wave inwave
	variable segment
		
	string wavestring = nameofwave(inwave)
	variable i,EndOfProtocol
	for (i = 1; i <= segment; i += 1)
		EndOfProtocol += GetSegmentDuration(wavestring,i)
	endfor
	
	return EndOfProtocol
end
//***********************************************************
static function EvaluateCalciumCurrents(ListOfWaves,RelevantSegment)
	string ListOfWaves
	variable RelevantSegment
	if(RelevantSegment < 0)
		DoAlert 0, "Non-existing segment"
		return -1
	endif
	string I_strings = ListMatch(ListOfWaves,"*Imon*")
	variable NumberOfWaves = Itemsinlist(I_strings)
	
	if(NumberOfWaves > 0)
		Make/O/T/N=(NumberOfWaves) I_traces = StringFromList(p,I_strings)
	endif

	variable StartOfProtocol,EndOfProtocol
	variable i
	for (i = 0; i < RelevantSegment; i += 1)
		StartOfProtocol += GetSegmentDuration(I_traces[1],i + 1)
	endfor
	EndOfProtocol = StartOfProtocol + GetSegmentDuration(I_traces[1],RelevantSegment + 1)

	Make/D/O/N=(NumberOfWaves) voltage,baseline, Ca_ss,Ca_Tail,Ca_activation	
	variable base, sdev, SteadyStateLevel
	
	for(i = 0; i < NumberOfWaves; i += 1)
		Duplicate/FREE $(I_traces[i]) WaveToAnalyse

		voltage[i] = str2num(stringbykey("V2",note($(I_traces[i])),":","\r"))
		wavestats /Q/R=(k_spacer,(StartOfProtocol - k_spacer))/M=2 WaveToAnalyse
		base = V_avg
		sdev = V_sdev
		baseline[i] = base * 1e12

		SteadyStateLevel = mean(WaveToAnalyse,(EndOfProtocol - k_spacer),EndOfProtocol )
		Ca_ss[i] = (SteadyStateLevel-base) * 1e12
		
		wavestats/Q/R=(EndOfProtocol, EndOfProtocol + k_spacer * 2) WaveToAnalyse
		Ca_Tail[i] = V_min * 1e12
		
		if (SteadyStateLevel < (base - 1 * sdev)) 
			variable V_fiterror
			CurveFit/NTHR=0/Q/N/W=2  exp_XOffset WaveToAnalyse[x2pnt(WaveToAnalyse,StartOfProtocol + k_spacer / 2),x2pnt(WaveToAnalyse,StartOfProtocol + k_spacer * 4)] /D 
			Ca_activation[i] = K2*1000
		else
			Ca_activation[i] = NaN
		endif
			
	endfor
	
	Duplicate/O Ca_Tail Ca_Tail_norm
	wavestats/Q Ca_Tail_norm
	Ca_Tail_norm /= V_min
	
	Edit/W=(605,25,1220,325) voltage,baseline, Ca_ss,Ca_Tail,Ca_Tail_norm,Ca_activation
	
	if (CheckIfWaveIsDisplayed(Ca_ss) == FALSE)
		Display /K=1 Ca_ss vs voltage
		Display /K=1 Ca_Tail_norm vs voltage
		Display /K=1 Ca_activation vs voltage
	
		Display /K=1
		for(i = 0; i < NumberOfWaves; i += 1)
			AppendtoGraph $stringfromlist(i,I_strings)
		endfor

		string cmd = "TileWindows/A=(3,2)/O=1/W=(0,0,600,720)"
		execute cmd
	endif
end

//***********************************************************
static function EvaluateDeactKinetics(ListOfWaves,RelevantSegment)
	string ListOfWaves
	variable RelevantSegment
	if(RelevantSegment < 0)
		DoAlert 0, "Non-existing segment"
		return -1
	endif
	string I_strings = ListMatch(ListOfWaves,"*Imon*")
	variable NumberOfWaves = Itemsinlist(I_strings)
	
	if(NumberOfWaves > 0)
		Make/O/T/N=(NumberOfWaves) I_traces = StringFromList(p,I_strings)
	endif


	variable StartOfProtocol,EndOfProtocol
	variable i
	for (i = 0; i < RelevantSegment; i += 1)
		StartOfProtocol += GetSegmentDuration(I_traces[1],i + 1)
	endfor
	EndOfProtocol = StartOfProtocol + GetSegmentDuration(I_traces[1],RelevantSegment + 1)

	Make/O/N=(NumberOfWaves) voltage,baseline,Ca_tail,Ca_deactivation	
	variable base, sdev, SteadyStateLevel	
	for(i = 0; i < NumberOfWaves; i += 1)
		Duplicate/FREE $(I_traces[i]) WaveToAnalyse

		voltage[i] = str2num(stringbykey("V3",note($(I_traces[i])),":","\r"))
		wavestats /Q/R=(k_spacer,(StartOfProtocol - k_spacer))/M=2 WaveToAnalyse
		base = V_avg
		sdev = V_sdev
		baseline[i] = base * 1e12
		
		wavestats/Q/R=(StartOfProtocol, StartOfProtocol + k_spacer * 2) WaveToAnalyse
		Ca_Tail[i] = V_min * 1e12
		
		variable V_fiterror
		CurveFit/NTHR=0/Q/N/W=2  exp_XOffset WaveToAnalyse[x2pnt(WaveToAnalyse,StartOfProtocol + k_spacer / 2),x2pnt(WaveToAnalyse,StartOfProtocol + k_spacer * 4)] /D 
		Ca_deactivation[i] = ((K2 > 0) && (K2 < 0.05)) ? K2*1000 : NaN		
	endfor
	
	Duplicate/O Ca_Tail Ca_Tail_norm
	wavestats/Q Ca_Tail_norm
	Ca_Tail_norm /= V_min
	
	Edit/W=(605,25,1220,325) voltage,baseline,Ca_Tail,Ca_deactivation
	
	if (CheckIfWaveIsDisplayed(Ca_deactivation) == FALSE)
		Display /K=1 Ca_deactivation vs voltage
	
		Display /K=1
		for(i = 0; i < NumberOfWaves; i += 1)
			AppendtoGraph $stringfromlist(i,I_strings)
		endfor

		string cmd = "TileWindows/A=(3,2)/O=1/W=(0,0,600,720)"
		execute cmd
	endif
	string results
	wfprintf results, "", Ca_deactivation
	putscraptext results
end

//***********************************************************
static function SubtractTraces()
	string prefix = "PMPulse_1_"
	string series1Str = "1"
	string series2Str = "2"
	string KillSourceWaves
	string analyse
	Prompt series1Str, "series to subtract from (" + prefix + "X_1_1):"
	Prompt series2Str, "series to subtract (" + prefix + "X_1_1):"
	Prompt KillSourceWaves, "Kill Source Waves",popup "yes;no"
	Prompt analyse, "Analyse subtracted traces",popup "yes;no"

	DoPrompt "Waves for Subtraction" series1Str,series2str,KillSourceWaves,analyse
	if (V_Flag)
		return -1
	endif

	KillLeakWaves()
	
	string ListOfWaves1 = WaveList(prefix+series1Str+"*",";","")
	if (itemsinlist(ListOfWaves1))
		wave TheWave1 = MakeMatrixWave(SelectImonTracesFromList(ListOfWaves1))
		Duplicate/O TheWave1 FirstMatrix
	else
		DoAlert 0, "No traces found."
		return -1
	endif

	string ListOfWaves2 = WaveList(prefix + series2Str + "*",";","")
	if (itemsinlist(ListOfWaves2))
		wave TheWave2 = MakeMatrixWave(SelectImonTracesFromList(ListOfWaves2))
		duplicate/O TheWave2 SecondMatrix
	else
		DoAlert 0, "No traces found."
		return -1
	endif

	if (DimSize(FirstMatrix,1) != dimsize(SecondMatrix,1) )
		DoAlert 0, "No. of waves mismatch."
		return -1
	endif

	Duplicate/O/FREE $(stringfromlist(0,ListOfWaves1)) DummyWave
	string OriginalWaveNote = note(DummyWave)
	MatrixOP/O SubtractedWaves = FirstMatrix - SecondMatrix
	Display/N=SubtractedTraces
	
	variable i
	variable NumOfWaves = dimsize(SubtractedWaves,1)
	for (i = 0; i < NumOfWaves; i += 1)
		AppendToGraph  SubtractedWaves[][i]
	endfor

	Killwaves/Z MatrixWave,FirstMatrix,SecondMatrix
	if(cmpstr(KillSourceWaves ,"yes") == 0)
		string WavesToKill = ListOfWaves1 + ";" + ListOfWaves2
		NumOfWaves = itemsinlist(WavesToKill)
		for(i = 0; i < NumOfWaves; i += 1)
			wave KillThisWave = $(stringfromlist(i,WavesToKill))
			KillWaves/Z KillThisWave
		endfor
	endif

	if(cmpstr(analyse ,"yes") == 0)
		NumOfWaves = dimsize(SubtractedWaves,1)
		for (i = 0; i < NumOfWaves; i += 1)
			wave ThisWave = $("K_subtracted" + "Imon" + num2str(i))
			MatrixOP/O ThisWave = col(SubtractedWaves,i)
			Note/K ThisWave, OriginalWaveNote
			CopyScales DummyWave ThisWave
		endfor
		
		KillWindow SubtractedTraces
		Killwaves/Z SubtractedWaves
		string ListOfWaves = WaveList("K_subtracted" + "*",";","")
		EvaluateCurrents(ListOfWaves,2,2)

		Display/N=SubtractedTraces
		
		variable NumOfItems = itemsinlist(ListOfWaves)
		for (i = 0; i < NumOfItems; i += 1)
			AppendToGraph $(stringfromlist(i,ListOfWaves))
		endfor
	endif

End

//***********************************************************
static function EvaluateTail(ListOfWaves,RelevantSegment)
	string ListOfWaves
	variable RelevantSegment
	if(RelevantSegment < 0)
		DoAlert 0, "Non-existing segment"
		return -1
	endif
	string I_strings = ListMatch(ListOfWaves,"*Imon*")
	variable NumberOfWaves = Itemsinlist(I_strings)
	
	if(NumberOfWaves > 0)
		Make/O/T/N=(NumberOfWaves) I_traces = StringFromList(p,I_strings)
	endif

	variable StartOfProtocol,EndOfProtocol
	variable i
	for (i = 0; i < RelevantSegment; i += 1)
		StartOfProtocol += GetSegmentDuration(I_traces[1],i + 1)
	endfor
	EndOfProtocol = StartOfProtocol + GetSegmentDuration(I_traces[1],RelevantSegment + 1)

	Make/O/N=(NumberOfWaves) TailCurrentAmp,TailCurrentIntegral
	
	variable base, NegativeLevel,PeakVal
	variable MakeVoltageWave = FALSE
	wave/D voltage
	if(waveexists(voltage) == FALSE)
		Make/O/D/N=(NumberOfWaves) voltage 
		MakeVoltageWave = TRUE
	endif	
	
	for(i = 0; i < NumberOfWaves; i += 1)
		Duplicate/FREE $(I_traces[i]) WaveToAnalyse
		
		if (MakeVoltageWave == TRUE)
			voltage[i] = str2num(stringbykey("V2",note($(I_traces[i])),":","\r"))
		endif
			
		base = mean(WaveToAnalyse,k_spacer,(StartOfProtocol - k_spacer))
		multithread WaveToAnalyse -= base
		
		NegativeLevel = -3e-11
		FindPeak/B=50/M=(NegativeLevel)/N/Q/R=(EndOfProtocol,EndOfProtocol + 0.002) WaveToAnalyse
		TailCurrentAmp[i] =  V_flag ? 0 : V_PeakVal * 1e12
		TailCurrentIntegral[i] = area(WaveToAnalyse,EndOfProtocol + 0.0001,EndOfProtocol + 0.001) * 1e12
	endfor
	wavestats/Q/M=1 TailCurrentAmp
	Duplicate/O TailCurrentAmp NormTailCurrent
	NormTailCurrent /= V_min
	
	if (WindowExists("results") == FALSE)
		Edit/N=results/W=(0,475,600,675) NormTailCurrent
	else
		AppendToTable/W=results NormTailCurrent
	endif	
	if (CheckIfWaveIsDisplayed(NormTailCurrent) == FALSE)
		Display/K=1 NormTailCurrent vs voltage
	endif

	string cmd = "TileWindows/A=(2,3)/O=1/W=(0,0,900,470)"
	execute cmd
	
	string clipboard
	wfprintf clipboard, "%g\t", NormTailCurrent
	putscraptext clipboard
end

//***********************************************************
static function/S SelectImonTracesFromList(ListOfWaveNames)
	string ListOfWaveNames
	string I_strings = ListMatch(ListOfWaveNames,"*Imon*")
	return I_strings
end

//***********************************************************
static function/Wave MakeMatrixWave(InputList)
	string InputList

	variable NofWaves = itemsinlist(InputList)
	string awave = stringfromlist(0,InputList)
	Make/FREE/D/O/N=(numpnts($awave),NofWaves) MatrixWave
		if(NofWaves > 1)
			Concatenate/O InputList,MatrixWave
		elseif(NofWaves == 1)
			wave tempwave = $(stringfromlist(0,InputList))
			MatrixWave = tempwave
		endif
	return MatrixWave
end

//***********************************************************
static function EvalPeak(InputWave,SegmentStart)
	wave InputWave
	variable SegmentStart

	variable PeakPos
	try
		PeakPos = FindPeakPosition(InputWave,SegmentStart,750,5)
	catch
		if( V_AbortCode == 100 )
			try
				V_AbortCode = 0
				PeakPos = FindPeakPosition(InputWave,SegmentStart,1500,5)
			catch
				if( V_AbortCode == 100 )
					try
						V_AbortCode = 0
						PeakPos = FindPeakPosition(InputWave,SegmentStart,1500,4)
					catch
						if( V_AbortCode == 100 )
							return 0
						endif
					endtry
				endif
			endtry
		endif
	endtry

	return InputWave[PeakPos]
end

//***********************************************************
static function FindPeakPosition(InputWave,SegmentStart,SmoothFactor,Threshold)
	wave InputWave
	variable SegmentStart,Threshold,SmoothFactor

	Duplicate/O InputWave SmoothedInput
	Smooth SmoothFactor, SmoothedInput
	Differentiate SmoothedInput/D=SmoothedInput_DIF
	Differentiate SmoothedInput_DIF/D=SmoothedInput_2ndDIF
	Wavestats/Q/R=(0,SegmentStart)SmoothedInput_2ndDIF
	variable SecDerThresh = Threshold * V_sdev

	variable SegmentEnd = SegmentStart + 0.015
	SegmentStart += 0.00025

	FindPeak/B=50/M=(SecDerThresh)/Q/P/R=(SegmentEnd,SegmentStart) SmoothedInput_2ndDIF
	AbortOnValue V_flag,100
	variable PeakPosition
	if(SmoothedInput_DIF[V_PeakLoc] < 0)
		FindLevel/EDGE=1/Q/P/R=[V_PeakLoc,V_PeakLoc + 200] SmoothedInput_DIF,0
		PeakPosition = V_LevelX
	elseif(SmoothedInput_DIF[V_PeakLoc] > 0)
		FindLevel/EDGE=1/Q/P/R=[V_PeakLoc,V_PeakLoc - 200] SmoothedInput_DIF,0
		PeakPosition = V_LevelX
	else
		PeakPosition = V_PeakLoc
	endif
	KillWaves/Z SmoothedInput,SmoothedInput_DIF,SmoothedInput_2ndDIF

	return PeakPosition
end

//***********************************************************
static function MakeRsCorrection(InputList)
	string InputList

	variable DoCorrection = 1
	variable CurrentType = 1

	Prompt DoCorrection, "Offline Rs correction:",popup "no;yes"
	Prompt CurrentType, "Ionic current:",popup "Na;K"

	DoPrompt "Rs correction" DoCorrection,CurrentType
	if (V_Flag)
		return -1
	endif
	if (DoCorrection == 1)
		return -1
	endif

	InputList = ListMatch(InputList,"*Imon*")
	string wavenote = note($(stringfromlist(0,InputList)))
	variable Rs,RSfraction,Cm
	sscanf stringbykey("RSeries",wavenote,":","\r"), "%f" , Rs
	sscanf stringbykey("rsFraction",wavenote,":","\r"), "%f" , RSfraction
	sscanf stringbykey("CSlow",wavenote,":","\r"), "%f" , Cm

	Rs -= Rs*RSfraction/100
	Rs *= 1e6
	Cm *= 1e-12

	variable reversal
	switch (CurrentType)
	case 1:
		reversal = 0.065
		break
	case 2:
		reversal = -0.105
	endswitch

	string ProtocolwaveList = replacestring("Imon1",InputList,"Stim-1")
	variable i
	variable NumOfWaves = itemsinlist(InputList)
	for(i = 0; i < NumOfWaves; i += 1)
		wave ThisDataWave = $(stringfromlist(i,InputList))
		wave ThisProtocolWave = $(stringfromlist(i,ProtocolwaveList))
		RsCorrection(ThisDataWave,ThisProtocolWave,Rs,Cm,1,1,reversal)
		Killwaves/Z ThisProtocolWave
	endfor
end
//***********************************************************
static function/S GetMetaData(wavenamestring)
	string wavenamestring
	
	String WaveMetadata = note($(wavenamestring))
	string filename
	filename = Parsefilepath(0,stringbykey("Filename",WaveMetaData,":", "\r"),"/", 1, 0)
	filename = removeending(filename)
	filename = removeending(filename, ".dat")
	
	variable CSlow = numberbykey("CSlow",WaveMetaData,":", "\r")
	variable Rseries = numberbykey("RSeries",WaveMetaData,":", "\r")
	variable Rscomp = numberbykey("rsFraction",WaveMetaData,":", "\r")
	variable cFast = numberbykey("cFastAmp1",WaveMetaData,":", "\r") + numberbykey("cFastAmp2",WaveMetaData,":", "\r")
	string Parameterstring
	sprintf Parameterstring, "\r%s\r%g\t%g\t%g\t%g",filename,Rseries,cFast,CSlow,Rscomp
	
	return Parameterstring
end	
//**********************************************************
static function GetSeriesNum()

	string AllWaves = WaveList("PM*",";","MINROWS:150")
	String TheWaveString = stringfromlist(0,AllWaves)
	
	string Prefix,Suffix
	variable v1,v2,v3,v4
	sscanf TheWaveString, "%[^_]%*[_]%d%*[_]%d%*[_]%d%*[_]%d%*[_]%s",prefix,v1,v2,v3,v4,Suffix
	
	return v2
end