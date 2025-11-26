; Copyright Pietro Terna, Gianpiero Pescarmona, Alberto Acquadro, Paolo Pescarmona, Giuseppe Russo, Emilio Sulis, Stefano Terna
; See Info tab for full copyright and license. Starting March 2020.
; Corresponding author pietro.terna@unito.it
;
; How to cite
;
; Terna P., Pescarmona G., Acquadro A., Russo G., Pescarmona P., Terna S. (2020), An Agent-Based Model of the Diffusione of Covid-19
; Using NetLogo, URL https://terna.to.it/simul/SIsaR.html
;
; Model website https://terna.to.it/simul/SIsaR.html
;
; A short paper at https://rofasss.org/2020/10/20/sisar/
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;WARNING

;The program contains the use of "int", "sort", "repeat", and "foreach" instructions, apparently useless.
;
;They are used to make the results consistent between the desk version and the web version.
;See my message of May 13, 2020, published at https://gitter.im/NetLogo/NetLogo,
;about three problems with the web version of NetLogo.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

__includes [ "vaccinatePeople.nls" ] ;comment for web version
extensions [csv] ;comment for web version

globals [csvList contagion# saveCSVfile infectingAgent
         setupPhase hLimit lLimit x aSchool aFactory rr pa pos chosen
         cumulativeRedCount         cumulativeRedCount_t-1
         cumulativeVioletCount      cumulativeVioletCount_t-1
         cumulativeGreenCount       cumulativeGreenCount_t-1
         cumulativeGreenCount*      cumulativeGreenCount_t-1*
         ; cumulativeGreenCount* and newGreenCount* only consider the symptomatic cases
         cumulativeDeadCount

         cumulativeRedCountNH
         cumulativeVioletCountNH
         cumulativeGreenCountNH
         cumulativeGreenCount*NH
         cumulativeDeadCountNH

         cumulativeRedCountStudents
         cumulativeVioletCountStudents
         cumulativeRedCountTeachers
         cumulativeVioletCountTeachers

        newRedCount newVioletCount newGreenCount newGreenCount* newDeadCount
         cityFrequency townFrequency houseFrequency candidatePatches countC countT
         countRepetitions limitations/Lockdown offPatches openPatch new-xcor new-ycor myBuilding radiusVector
         numberNursingHomes numberOperatorsPerNH probOfBeingInNursingHome% numberHospitals visitingNHsP%
         mySchool schoolsOn freePatches TinfSch SinfSch TinfNotSch SinfNotSch studentNumber teacherNumber
         myFactory factoriesOn largeUnits smallUnits workerNumberRegular workerNumberFragile quotaRegularWorkers quotaFragileWorkers
         open open-1 tot isStud?

         checkPoint1activeAt checkPoint2activeAt checkPoint3activeAt checkPoint4activeAt checkPoint5activeAt checkPoint6activeAt
;        these are too long nema, but explicative, so we keep them as comments, using those below to have narrow columns in Behavior Space tables
;        with chekPoint1cumulativeRedCount cp1cRed etc.
;         chekPoint1cumulativeRedCount chekPoint1cumulativeRedCountNH chekPoint1cumulativeVioletCount chekPoint1cumulativeVioletCountNH chekPoint1cumulativeDeadCount
;         chekPoint2cumulativeRedCount chekPoint2cumulativeRedCountNH chekPoint2cumulativeVioletCount chekPoint2cumulativeVioletCountNH chekPoint2cumulativeDeadCount
;         chekPoint3cumulativeRedCount chekPoint3cumulativeRedCountNH chekPoint3cumulativeVioletCount chekPoint3cumulativeVioletCountNH chekPoint3cumulativeDeadCount
         cp1cRed cp1cRedNH cp1cViolet cp1cVioletNH cp1cDead
         cp2cRed cp2cRedNH cp2cViolet cp2cVioletNH cp2cDead
         cp3cRed cp3cRedNH cp3cViolet cp3cVioletNH cp3cDead
         cp4cRed cp4cRedNH cp4cViolet cp4cVioletNH cp4cDead
         cp5cRed cp5cRedNH cp5cViolet cp5cVioletNH cp5cDead
         cp6cRed cp6cRedNH cp6cViolet cp6cVioletNH cp6cDead

         myScript lll howManyRows s n l nsub sl i choice seed seed0 seedD; in this row, tools for script decode and use
         displayStatus howManyDecorativeAgents peopleSize0 z tmpAg
         theEnd
         delayedStartList countD
         iw fileVacc useGA
         fileName outList0
        ]

breed [people person]
breed [decorativeAgents aDecAg]

people-own [myContagion# whereInfectedColor myFragilityRate myInfectingAgent myInfectingAgentContagion# myInfectingAgentWhereInfectedColor contagionTick
            symptomaticAsymptomaticStatus
            startingInfection finishingInfection deadTick xcor-1 ycor-1 xcorHome ycorHome apartment
            delayedStart delayedFinish
            student teacher xcorSchool ycorSchool classroom ; xcorSchool ycorSchool classroom set to -1 if not student teacher OR the school is not active
            inNursingHome infectedInHospital infectedAtHome intrinsicSusceptibility young temporaryDestination
            subjectiveRandomValueInRepetitions #visitsInNHs infectedAtSchool infectedAtWork worker xcorFactory ycorFactory vaccinationDate activeVaccination
            contagiousAfterVaccineActivation myGroupForStatistics]

patches-own [sdReductionF off downtown dist teacherCount studentCount workerCount openFactory]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to setup
  clear-all
  reset-ticks
  set setupPhase True

  decodeScript
  checkScript

  setPars&Vars

  createWorld
  createPeople


  createInitialInfected

  createDelayedInfected

  if fileVacc != 0 [vaccinatePeople set outList0 outSusceptibleList0] ;comment for web version

  set setupPhase False
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go
  tick
  clearData

  checkScript
  makeCheckPoints

  checklimitations/Lockdown
  checkBufferZones

  checkPeopleSize
  setDisplay

  activateDelayedInfections

  sendPeopleHome
  diffuseInfectionAtHomeOrInNursingHomesOrInHospitalsNightStep

  set countRepetitions 0
     while [countRepetitions < #movementCyclesPerTick]
     [

      makeMovement

     if countRepetitions = 0
      [

       if fileVacc != 0 [haveVaccinatedPeopleComingOut] ;comment for web version
       evolveInfectedPeople    ; to Recovered
       haveInfectedPeopleComingOut ; infected in previous periods and that have finished incubation
       diffuseInfectionWithinFactories
       diffuseInfectionAtSchool

      ]

      diffuseInfectionExcludingFactoriesSchools

      set countRepetitions countRepetitions + 1
     ]

  if fileVacc != 0 [cleanInfections] ;comment for web version

  makePlots

  ;to finish

    if
     useGA = "n" and
     count people with [startingInfection + 1 > ticks]  = 0 and
     count people with [finishingInfection + 1 > ticks] = 0 and
     count people with [delayedStart > 0] = 0                  [set theEnd True] ; used also in BehaviorSpace

    if theEnd
     [makeFinalOutputAndCalculations
      stop]

  if myStop != 0 and ticks = myStop [stop]

  if fileVacc != 0 and ticks = vaccinationDay [set susceptibleListV countSusceptible] ;comment for web version

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to setPars&Vars

  ; controlling with GA
  set useGA "n"

  ; concluding, here and in BehaviorSpace
  set theEnd False

  ; for the csv file of contagions
  set csvList []

  ; contagion progressive number
  set contagion# 1

  ; schools still missing
  set schoolsOn False

  ; factories still missing
  set factoriesOn False
  set largeUnits 5    ; look at Info sheet
  set smallUnits 100  ; look at Info sheet

  ; workers, look at Info view
  set quotaRegularWorkers 0.86
  set quotaFragileWorkers 0.16

  ;limits of the space
  set hLimit max-pycor
  if not cities&towns [set hLimit 50] ; hospitals & nunsing homes  bottom limit 50

  set lLimit min-pycor
  if not cities&towns [set lLimit 3] ; house zone top limit 3 (4 rows from 0 to 3)

  ; how much far from home
  set radiusVector [2 3 3 3 3 5 5 5 7 9]

  ;cities and towns, if any
  set cityFrequency 0.0025
  set townFrequency 0.0125
  set houseFrequency 0.5

  set numberNursingHomes 5
  set numberOperatorsPerNH 3
  set probOfBeingInNursingHome% 5 ; as a quota of fragile people, explanation in Info, sub nursing homes
  set visitingNHsP% 1.5 ; look at Info

  set numberHospitals 5

  ; zeroing
  set cumulativeRedCount     0 set cumulativeVioletCount     0 set cumulativeGreenCount     0
  set cumulativeRedCountNH   0 set cumulativeVioletCountNH   0 set cumulativeGreenCountNH   0
  set cumulativeRedCount_t-1 0 set cumulativeVioletCount_t-1 0 set cumulativeGreenCount_t-1 0
  set cumulativeDeadCount    0 set cumulativeDeadCountNH     0
  set cumulativeRedCountStudents 0 set cumulativeVioletCountStudents 0
  set cumulativeRedCountTeachers 0 set cumulativeVioletCountTeachers 0

  ;the gray color is represented by gray + 2, lighter than standard gray
  ;the green color is set to green + 1
  ;the turquoise color is set to turquoise - 1
  ;the pink color is represented by pink + 3, lighter than standard pink for hospitals not activated
  ;but the active hospitals are pink
  ;the cyan color is represented by cyan + 3, lighter than standard cyan
  ;the orange color is represented by orange + 3, lighter than standard orange for nursing homes not activated
  ;but the active nurging homes are orange
  ;into the plots is represented we use the lime color
  ;the brown color is used as brown - 2 for large production units and brown + 2 for small ones

  ; check points for Behavior Space

  set checkPoint1activeAt -1
  set checkPoint2activeAt -1
  set checkPoint3activeAt -1
  set checkPoint4activeAt -1
  set checkPoint5activeAt -1
  set checkPoint6activeAt -1
  set cp1cRed -1 set cp1cRedNH -1 set cp1cViolet -1 set cp1cVioletNH -1 set cp1cDead -1
  set cp2cRed -1 set cp2cRedNH -1 set cp2cViolet -1 set cp2cVioletNH -1 set cp2cDead -1
  set cp3cRed -1 set cp3cRedNH -1 set cp3cViolet -1 set cp3cVioletNH -1 set cp3cDead -1
  set cp4cRed -1 set cp4cRedNH -1 set cp4cViolet -1 set cp4cVioletNH -1 set cp4cDead -1
  set cp5cRed -1 set cp5cRedNH -1 set cp5cViolet -1 set cp5cVioletNH -1 set cp5cDead -1
  set cp6cRed -1 set cp6cRedNH -1 set cp6cViolet -1 set cp6cVioletNH -1 set cp6cDead -1

  set limitations/Lockdown False

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to createWorld

  ask patches [set pcolor white set off False set downtown False
               set teacherCount 0 set studentCount 0 set workerCount 0 set openFactory False]
  ; setting movement standard deviation into the patches
  ask patches [set sdReductionF 1]

  ;---------

  ;hospital or nursingHomes, potential area
  if not cities&towns
    [ask patches with [pycor >= hLimit and pxcor <= 23][set pcolor pink + 3]
     ask patches with [pycor >= hLimit and pxcor >= 25][set pcolor orange + 3]

      if activateNursingHomes [ask n-of numberNursingHomes patches with [pcolor = orange + 3]
         [set pcolor orange]] ; create active nursing homes

      if activateHospitals    [ask n-of numberHospitals    patches with [pcolor = pink + 3]
         [set pcolor pink]]] ; create active hospitals

  ;separation
  if not cities&towns
    [ask patches with [pycor >= hLimit and pxcor = 24][set pcolor brown]]

  ;houses
  if not cities&towns
    [ask patches with [pycor <= lLimit][set pcolor cyan + 3]]  ; about 8% of patches with houses, accordin to ($) below

  ;---------

  if cities&towns [
    set countC 0 set countT 0
    while [countC < 3 or countT < 6]
    [
       ask patches [let f random-float 1
                    if f < cityFrequency [set plabel "C" set plabel-color blue - 2 set downtown True]
                    if f >= cityFrequency and f < (townFrequency + cityFrequency)
                                         [set plabel "T" set plabel-color blue - 2 set downtown True]
                   ]
      ask patches with [plabel = "T"]   [ask patches in-radius 1.5 [set plabel "T" set plabel-color blue - 2]]
      ask patches with [plabel = "C"]   [ask patches in-radius 3   [set plabel "C" set plabel-color blue - 2]]

      set countC count patches with [plabel = "C"]
      set countT count patches with [plabel = "T"]
    ]
  ]

  if cities&towns [
     ask patches with [plabel = "T"][set sdReductionF 2]
     ask patches with [plabel = "C"][set sdReductionF 3]
               ]


  ;houses
  if cities&towns  ; explantions in Info
    [ask patches with [(plabel = "C" or plabel = "T")]
                                                     [if random-float 1 < houseFrequency [set pcolor cyan + 3]]] ; 1/2 of the patches with C or T have houses


  ;schools
  if cities&towns and activateSchools [
    set schoolsOn True
    ask patches with [plabel = "T" and downtown][if random 5 >= 3  [set freePatches patch-set sort patches in-radius 1.5 with [pcolor = white]
                                                                    if any? freePatches
                                                                    [set aSchool one-of sort freePatches; 2 schools for 5 towns
                                                                     if aSchool != nobody [ask aSchool [set pcolor yellow]]]]]
    ask patches with [plabel = "C" and downtown][if random 2 >= 0  [set freePatches patch-set sort patches in-radius 2    with [pcolor = white]  ; 1 schools for each city
                                                                    if any? freePatches
                                                                    [set aSchool one-of sort freePatches; 1 schools for each city
                                                                     if aSchool != nobody [ask aSchool [set pcolor yellow]]]]];random is kept for other possible calibrations
    ; sort above are ncessary to mantain consinstency betwewn the desk and the web versions of the program
  ]


  ;hospital or nursingHomes

  if cities&towns and activateNursingHomes
   [set candidatePatches patches with [pcolor = white and not (plabel = "C" or plabel = "T")]              ; other candidates
   ; numberNursingHomes, as we explain in Info
      ask n-of numberNursingHomes candidatePatches with [count patches in-radius 5 with [downtown] > 0] ; not so far from urban areas
      [set pcolor orange set plabel ""]]


  if cities&towns and activateHospitals
    [set candidatePatches patches with [pcolor = white and (plabel = "C" or plabel = "T")]         ; candidate patches
    ask n-of numberHospitals candidatePatches [set pcolor pink set plabel ""]]                     ; numberHospitals, as we explain in Info


  if cities&towns [
     ask patches with [pcolor = orange] [set sdReductionF 2] ; coming from no label patches, so the same value of "T"
     ask patches with [pcolor = pink  ] [set sdReductionF 2] ; set to a very low value, independerly from the origin "C" o "T"


               ]
  ;factories
  if cities&towns and activateFactories
   [
    set factoriesOn True

    set freePatches patch-set n-of largeUnits sort patches with [(plabel = "T" or plabel = "C")]
    ask freePatches [set tot 0 foreach (sort patches in-radius 10 with [pcolor = white and plabel = ""
                                                     and count neighbors with [pcolor = orange or pcolor = pink or pcolor = yellow] = 0])[set tot tot + 1]
                     if tot > 0 [set aFactory one-of sort patches in-radius 10 with [pcolor = white and plabel = ""
                                              and count neighbors with [pcolor = orange or pcolor = pink or pcolor = yellow] = 0] ; for largeUnits value, look at the Info sheet
                     ;if aFactory != nobody [ask aFactory [set pcolor brown - 2]]]; the nobody testi is unuseful, never nobody because one-of gives nobody on an empty agentset and error on an empy list
                     ask aFactory [set pcolor brown - 2]]]
                     ; large units

    set freePatches patch-set n-of smallUnits sort patches with [(plabel = "T" or plabel = "C")]
    ask freePatches [set tot 0 foreach (sort patches in-radius 5 with [pcolor = white and plabel = ""
                                                    and count neighbors with [pcolor = orange or pcolor = pink or pcolor = yellow] = 0]) [set tot tot + 1]
                     if tot > 0 [set aFactory one-of sort patches in-radius 5  with [pcolor = white and plabel = ""
                                                and count neighbors with [pcolor = orange or pcolor = pink or pcolor = yellow] = 0] ; for largeUnits value, look at the Info sheet
                     ;if aFactory != nobody [ask aFactory [set pcolor brown + 2]]]; the nobody testi is unuseful, never nobody because one-of gives nobody on an empty agentset and error on an empy list
                     ask aFactory [set pcolor brown + 2]]];
                     ; small units
   ]
    ; sort above are ncessary to mantain consinstency betwewn the desk and the web versions of the program



end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to createPeople
  create-people populationSize
         [
            set myContagion# 0
            set whereInfectedColor 0
            set myFragilityRate 0
            set myInfectingAgent 0
            set myInfectingAgentContagion# 0
            set myInfectingAgentWhereInfectedColor 0
            set contagionTick 0
            set symptomaticAsymptomaticStatus 0

            set color gray + 2
            set shape "cylinder"
            set size 0.4
            set startingInfection 0
            set finishingInfection 0
            set deadTick 0

            set delayedStart 0
            set delayedFinish 0

            assignPeopleHome
            movetoWhiteCyan ; everyone; people in nursingHome will have their home re-assigned below and will be there

            set xcor-1 0 set ycor-1 0
            set inNursingHome False
            set infectedInHospital False
            set infectedAtHome False
            set intrinsicSusceptibility 1 / intrinsicSusceptibilityFactor ; for regular people
            set myFragilityRate 2
            set young False
            set temporaryDestination False
            set subjectiveRandomValueInRepetitions 0
            set student False
            set teacher False
            set xcorSchool -1
            set ycorSchool -1
            set classroom -1
            set worker False
            set xcorFactory -1
            set ycorFactory -1
            set #visitsInNHs 0
            set infectedAtSchool False
            set infectedAtWork False
            set vaccinationDate 0
            set activeVaccination false
            set contagiousAfterVaccineActivation true
            set myGroupForStatistics 0

         ]

  ; fragile and young people
  ask people [if random-float 100 <= fragilePeople% [set shape "x" set size 0.8
                                     set intrinsicSusceptibility 1
                                     set myFragilityRate 3]] ; for fragile people
  ; on residuals, youngPeople% is x in population*(1-fragilePeople%)*x = population*youngPeople%
  set x youngPeople% / (1 - fragilePeople% / 100)
  ask people with [shape = "cylinder"][if random-float 100 <= x [set young True set myFragilityRate 1
                                                               set intrinsicSusceptibility 1 / intrinsicSusceptibilityFactor ^ 2]]

  ; healtcare operators (hospital related) in population not fragile
  if activateHoperators  [
  ask people with [shape = "cylinder"] [if random-float 100 <= healthcareOperators% [set shape "star" set size 1]]]

  ; healtcare operators (nursing home related) in population not fragile
  if activateNHoperators [
  let q numberOperatorsPerNH * numberNursingHomes
  ask n-of q people with [shape = "cylinder"] [set shape "triangle" set size 1]]


  ; people in nursing homes
      ask people with [shape = "x"] [if random-float 100 < probOfBeingInNursingHome%
                                    [let nNH count patches with [pcolor = orange]
                                     if  nNH > 0 [let myNH one-of patches with [pcolor = orange]
                                                  set xcor [pxcor] of myNH
                                                  set ycor [pycor] of myNH
                                                  set xcorHome xcor set ycorHome ycor set apartment 0
                                                  set xcor xcor  - 0.3 + random-float 0.6
                                                  set ycor ycor  - 0.3 + random-float 0.6
                                                  set xcor-1 xcor set ycor-1 ycor
                                                  set inNursingHome True
                                                  set intrinsicSusceptibility intrinsicSusceptibilityFactor
                                                  set myFragilityRate 4] ; for extra-fragile people

                                     ]
                                   ]

  ; theachers and students
       set teacherNumber 4 * count patches with [pcolor = yellow]
       ask n-of teacherNumber people with [shape = "cylinder" and not young][set teacher True]

       set studentNumber int (2 * count people with [young] / 3)
       ask n-of studentNumber people with [shape = "cylinder" and young][set student True]

  if cities&towns and activateSchools
      [
       ; teachers and closer school with classroom
       ask people with [teacher]
               [
                set xcor-1 xcor set ycor-1 ycor
                setxy xcorhome ycorhome ; temporary
                ask patches with [pcolor = yellow and teacherCount <= 3] ; so 4 is maximum
                       [set dist  distance myself]
                set myschool first sort-on [dist]  patches with [pcolor = yellow and teacherCount <= 3]
                  ask myschool [set teacherCount teacherCount + 1]
                set xcorSchool [pxcor] of myschool
                set ycorSchool [pycor] of myschool
                set classroom 1 + random howManyClassroomsInASchool ; can be 1 or 2
                setxy xcor-1 ycor-1
               ]

       ; students and closer school with classroom
       let limit int (count people with [student] / count patches with [pcolor = yellow]) + 2
       ask people with [student]
               [
                set xcor-1 xcor set ycor-1 ycor
                setxy xcorhome ycorhome ; temporary
                ask patches with [pcolor = yellow and studentCount <= limit]
                       [set dist  distance myself]
                set myschool first sort-on [dist]  patches with [pcolor = yellow and studentCount <= limit]
                  ask myschool [set studentCount studentCount + 1]
                set xcorSchool [pxcor] of myschool
                set ycorSchool [pycor] of myschool
                set classroom 1 + random howManyClassroomsInASchool ; can be 1 or 2
                setxy xcor-1 ycor-1
               ]
     ]


 ; workers, look in Info, workplaces and employees
  if not fragileWorkersAtHome
    [
     set workerNumberRegular int (quotaRegularWorkers * count people with [shape = "cylinder" and not young  and not teacher])
     set workerNumberFragile int (quotaFragileWorkers * count people with [shape = "x"])
     ask n-of workerNumberRegular people with [shape = "cylinder" and not young and not teacher][set worker True]
     ask n-of workerNumberFragile people with [shape = "x" and not inNursingHome][set worker True]
    ]
  if fragileWorkersAtHome ; no fragile, same total
    [
     set workerNumberRegular int (quotaRegularWorkers * count people with [shape = "cylinder" and not young  and not teacher])
     set workerNumberFragile int (quotaFragileWorkers * count people with [shape = "x"])
     ask n-of (workerNumberRegular + workerNumberFragile) people with [shape = "cylinder" and not young and not teacher][set worker True]
    ]

  if cities&towns and activateFactories
      [
       ; workers and closer factory
       ask people with [worker]
               [
                set xcor-1 xcor set ycor-1 ycor
                setxy xcorhome ycorhome ; temporary
                ifelse random-float 3 <= 1
                  [ask patches with [pcolor = brown - 2] ; large factories
                                    [set dist distance myself]
                   set myfactory first sort-on [dist]  patches with [pcolor = brown - 2 and workerCount < 150]
                   ask myfactory [set workerCount workerCount + 1]
                  ]
                  [ask patches with [pcolor = brown + 2] ; small factories
                                    [set dist distance myself]
                   set myfactory first sort-on [dist]  patches with [pcolor = brown + 2 and workerCount < 15]
                   ask myfactory [set workerCount workerCount + 1]
                  ]
                set xcorFactory [pxcor] of myfactory
                set ycorFactory [pycor] of myfactory
                setxy xcor-1 ycor-1
               ]
      ]



  ; to avoid people aligned in the center of the patches and to move people in nursing home in their n. h.
  repeat 2   [makeMovement]

  ; preparing 3D with decorative agents (turtles not moving, with special shapes)
  ;ask one-of patches [set pcolor brown] ; test
  ask patches with [pcolor = cyan + 3]     [sprout-decorativeAgents 1 [set shape "house"                     set color cyan + 3  set size 1.5]]
  ask patches with [pcolor = orange]       [sprout-decorativeAgents 1 [set shape "house colonial"            set color orange    set size 2.0]]
  ask patches with [pcolor = pink]         [sprout-decorativeAgents 1 [set shape "building institution"      set color pink      set size 2.0]]
  ask patches with [pcolor = yellow]       [sprout-decorativeAgents 1 [set shape "house ranch"               set color yellow    set size 2.0]]
  ask patches with [pcolor = brown - 2]    [sprout-decorativeAgents 1 [set shape "factory"                   set color brown - 2 set size 3.0]]
  ask patches with [pcolor = brown + 2]    [sprout-decorativeAgents 1 [set shape "factory"                   set color brown - 2 set size 2.0]]
  ask decorativeAgents [set hidden? True]

  set displayStatus "DisplayPeople"
  set displayChoice "DisplayPeople"

  set peopleSize0 "Regular"
  set peopleSize  "Regular"

  set howManyDecorativeAgents count decorativeAgents


end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to createInitialInfected
  ask n-of int (ratioInitialInfected% * populationSize / 100) people with [pcolor != pink and pcolor != orange] ; not in hospital or nursing home
                                       [
                                        set startingInfection  1
                                        set finishingInfection startingInfection + minInfectionDuration
                                            + random (maxInfectionDuration - minInfectionDuration)

                                        ;output seguence prepaparion (set below or ";" already set or ";;" to be set in haveInfectedPeopleCommingOut:
                                          set myInfectingAgent 0                   ;external unindified agent
                                          set myInfectingAgentContagion# 0         ;as above
                                          set myInfectingAgentWhereInfectedColor 0 ;as above
                                          ;who
                                          ;myFragilityRate
                                          set myContagion# contagion# set contagion# contagion# + 1
                                          set whereInfectedColor 0                 ;as above
                                          set contagionTick ticks
                                          ;startingInfection
                                          ;finishingInfection
                                          ;;symptomaticAsymptomaticStatus

                                        ]
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to createDelayedInfected

  if length delayedStartList > 0 [

  if seedD = 0 [user-message "To use add@ you have to define seedD with 0 seedD 'some value'" stop]

  ;show seedD
  if seedD = -9999 [set seedD int (timer * 1000)]
  if seedD != -9999 and seedD < 0 [set seedD fromSeedSequence abs seedD]
  ;show seedD

  random-seed seedD

  ;show delayedStartList

  ask n-of length delayedStartList people with [pcolor != pink and pcolor != orange and startingInfection = 0]
     ; not in hospital or nursing home and not in the firt infected set
     [;show who
      set delayedStart first delayedStartList + incubationPeriod
      set delayedStartList but-first delayedStartList
      set delayedFinish delayedStart + minInfectionDuration
                        + random (maxInfectionDuration - minInfectionDuration)
     ]
  ]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to checklimitations/Lockdown

  if activateLimitations/LockdownFromTick = ticks
      [
        set limitations/Lockdown True]

  if activateLimitations/LockdownFromTick = 0
      [
        set limitations/Lockdown False]


end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to hardFinishLockdown

  set activateLimitations/LockdownFromTick 0 checklimitations/Lockdown

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to sendPeopleHome

       ; people in nursing homes have their home there
       ask people
        [
          if not infectedInHospital and patch xcorHome ycorHome != patch-here [
                                    if pcolor != pink and not temporaryDestination ; back from a H or directly from a visit to a NH or
                                                                                   ; from school: keep old -1 positions
                                              [set xcor-1 xcor set ycor-1 ycor]
                                               setxy xcorHome - 0.3 + random-float 0.6 ycorHome  - 0.3 + random-float 0.6]
                                    set temporaryDestination False
                                    ; to switch off temporary status if coming directly from visiting NHs or going to school
       ]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to makeMovement

  if cities&towns and activateSchools and ticks > 0 and not schoolsOn
     [user-message "It is impossible to activate schools in a tick > 0 if the activation at tick 0 is missing." stop]

  ; check open factories
  if not setupPhase and activateFactories
  [
   set open int (0.01 * %openFactoriesAlsoWhenLimitationsOn * count patches with [pcolor = (brown - 2) or pcolor = (brown + 2)])
   set open-1 count patches with [openFactory]
   if open > open-1 [ask n-of (open - open-1) patches with [not openFactory and
                                                     (pcolor = (brown - 2) or pcolor = (brown + 2))][set openFactory True]]
   if open < open-1 [ask n-of (open-1 - open) patches with [openFactory][set openFactory False]]
  ]


  ask people      [if countRepetitions = 0 [set subjectiveRandomValueInRepetitions random-float 100]
                   ; to assure the consinstence among cnoices in the same tick, in "people can move if"

                   ; infected people sent to hospital or left at home (infectedInHospital or infectedAtHome)

                   ; case no hospitals
                   if color = red and not infectedInHospital and not infectedAtHome and not activateHospitals
                       [
                        set infectedAtHome True
                        ; already at home: if this is also the previous position, create a future destination to restart at the end of the infection
                        if patch xcorHome ycorHome = patch-here and patch xcor-1 ycor-1 = patch-here
                                                                          [let futureDestination one-of patches with [pcolor = white]
                                                                           set xcor-1 [pxcor] of futureDestination - 0.3 + random-float 0.6
                                                                           set ycor-1 [pycor] of futureDestination - 0.3 + random-float 0.6]
                       ; not at home, go
                       if patch xcorHome ycorHome != patch-here
                                                                          [set xcor-1 xcor
                                                                           set ycor-1 ycor
                                                                           set xcor xcorhome - 0.3 + random-float 0.6
                                                                           set ycor ycorhome - 0.3 + random-float 0.6
                                                                          ]
                      ]

                   ; case with hospitals
                   if  color = red and not infectedInHospital and not infectedAtHome and activateHospitals
                       [
                        ifelse random-float 100 < assignmentRateInfectedToHospitals%
                        ; infected sent in hospital
                          [
                          ;keep xcor-1 ycor-1 to be used after going back home as Recovered; if previously she/he was at home, create a future destination to restart at the end of the infection
                          if patch xcor-1 ycor-1 = patch xcorhome ycorhome [let futureDestination one-of patches with [pcolor = white]
                                                                           set xcor-1 [pxcor] of futureDestination - 0.3 + random-float 0.6
                                                                           set ycor-1 [pycor] of futureDestination - 0.3 + random-float 0.6]

                          let destination one-of patches with [pcolor = pink]
                          setxy [pxcor] of destination  - 0.3 + random-float 0.6
                                [pycor] of destination  - 0.3 + random-float 0.6
                          set infectedInHospital True
                         ]
                        ; infected sent at home
                         [
                          set infectedAtHome True
                          ; already at home: if this is also the previous position, create a future destination to restart at the end of the infection
                          if patch xcorHome ycorHome = patch-here and patch xcor-1 ycor-1 = patch-here
                                                                          [let futureDestination one-of patches with [pcolor = white]
                                                                           set xcor-1 [pxcor] of futureDestination - 0.3 + random-float 0.6
                                                                           set ycor-1 [pycor] of futureDestination - 0.3 + random-float 0.6]
                          ; not at home, go
                          if patch xcorHome ycorHome != patch-here
                                                                          [set xcor-1 xcor
                                                                           set ycor-1 ycor
                                                                           set xcor xcorhome - 0.3 + random-float 0.6
                                                                           set ycor ycorhome - 0.3 + random-float 0.6]
                        ]
                       ]

;******************************************************************************

                   ; people can move if

                                         ; RULE:  no limitations/Lockdown, not in hospital, not in nursing home, no symptomatic (controlled above)
                      if (not limitations/Lockdown and not inNursingHome and not infectedInHospital and not infectedAtHome) or

                                         ; EXCEPTIONS
                      (shape = "star" and color != red) or
                                         ; hospital healthcare operators not symptomatic

                      (shape = "triangle" and color != red) or
                                         ; nursing home healthcare operators not symptomatic
;
                      (shape != "x"   and color != red and not inNursingHome and subjectiveRandomValueInRepetitions < %PeopleNotFragileNotSymptomaticLeavingHome) or
                                         ; no fragile, no nursing home, no symptomatic with prob

                      (color != red and not inNursingHome and subjectiveRandomValueInRepetitions < %PeopleAnyTypeNotSymptomaticLeavingHome) or
                                         ; any type not symptomatic with prob

                      (color != red and worker and cities&towns and activateFactories and countRepetitions = 0 and
                                       [openFactory] of patch xcorFactory ycorFactory
                                       and not (stopFragileWorkers and shape = "x") ) or
                                         ; workers and open factories

                      (color != red and (student or teacher) and cities&towns and activateSchools and countRepetitions = 0)
                                         ; students and teachers and school active

;******************************************************************************

    [
                    set chosen False

;----------------------------------------
                    ; at home? if yes, go to the last "day before" place, where you were before going home (this is rule $, see below)
                    if patch xcorHome ycorHome = patch-here
                        [
                         setxy xcor-1 ycor-1
                         ]
;----------------------------------------

                    set new-xcor -1 set new-ycor -1

;----------------------------------------
                    ; go to visit a nursing home
                    if not chosen and peopleVisitingNHs and countRepetitions = 0 and activateNursingHomes and random-float 100 < visitingNHsP%                                  and not activateBufferZones
                                  and shape != "star" and shape != "triangle" and not setupPhase ;visiting a nursing home or making a regular movement
                       [
                        set chosen True
                        set #visitsInNHs #visitsInNHs + 1
                        let destination one-of patches with [pcolor = orange]
                        set new-xcor [pxcor] of destination  - 0.3 + random-float 0.6
                        set new-ycor [pycor] of destination  - 0.3 + random-float 0.6
                        set temporaryDestination True
                        setxy new-xcor new-ycor
                       ]

;----------------------------------------
                    ; go to a factory
                    if not chosen and worker and (cities&towns and activateFactories) and countRepetitions = 0
                                  and not setupPhase and not (stopFragileWorkers and shape = "x") ;going to work
                       [
                        if [openFactory] of patch xcorFactory ycorFactory
                        [
                         set chosen True
                         set temporaryDestination True
                         setxy xcorFactory - 0.3 + random-float 0.6 ycorFactory - 0.3 + random-float 0.6
                        ]
                       ]

;----------------------------------------
                    ; go to school
                    set isStud? true
                    if not student [set isStud? false]     ; to avoid changing old random sequences, choice is splitted
                    if student and %Students != 100 [if random-float 100 > %Students [set isStud? false]]
                    if not chosen and (isStud? or teacher)
                                  and (cities&towns and activateSchools) and countRepetitions = 0
                                  and not setupPhase ;going to school
                       [
                        set chosen True
                        set temporaryDestination True
                        setxy xcorSchool - 0.3 + random-float 0.6 ycorSchool - 0.3 + random-float 0.6
                       ]

;----------------------------------------
                    ; H operators going to work or moving around to any place
                    if not chosen and shape = "star"
                       [
                        set chosen True

                        ; work day or free day
                        if countRepetitions = 0
                        [
                        ifelse random-float 100 < 70; 5 working days per week, so around 70% e 30%
                          ; go to work (if not already there, due to rule ($) above)
                          [if activateHospitals and [pcolor] of patch-here != pink
                              [
                               let destination one-of patches with [pcolor = pink]
                               set new-xcor [pxcor] of destination  - 0.3 + random-float 0.6
                               set new-ycor [pycor] of destination  - 0.3 + random-float 0.6
                               set xcor-1 xcor set ycor-1 ycor
                               setxy new-xcor new-ycor
                              ]
                          ]
                          ; move around to any place, but id already in a white or cyan + 3 place make short random movements
                          [ifelse [pcolor] of patch-here = pink

                           [let destination one-of patches with [pcolor = white or pcolor = cyan + 3]
                            setxy [pxcor] of destination - 0.3 + random-float 0.6
                                  [pycor] of destination - 0.3 + random-float 0.6]

                           [set iw 0
                            while [new-xcor < (min-pxcor - 0.499) or new-xcor > (max-pxcor + 0.499) or
                               new-ycor < (min-pycor - 0.499) or new-ycor > (max-pycor + 0.499)
                               or [off] of patch new-xcor new-ycor = True or [pcolor] of patch new-xcor new-ycor = yellow
                              ]
                               [set iw iw + 1 if iw > 100000 and inspace and [pcolor] of patch new-xcor new-ycor = yellow [makeJump xcor ycor]
                                ; the control above is to avoid situation of perpetual while which, sometime, occurred here, but also to preserve
                                ; long repetition cycles effectively finishing and so  determining the random values sequence (in a back-compatibility effort)
                                if inSpace and [off] of patch-here = True [makeJump xcor ycor]
                                set new-xcor xcor + random-normal 0 sdRandomMovement / sdReductionF set new-ycor ycor + random-normal 0 sdRandomMovement / sdReductionF]

                           set xcor-1 xcor set ycor-1 ycor
                           setxy new-xcor new-ycor]
                           ]
                          ]
                        ]


;----------------------------------------
                    ; NH operators going to work or moving around to any place
                    if not chosen and shape = "triangle"
                       [
                        set chosen True

                        ; work day or free day
                        if countRepetitions = 0
                        [
                        ifelse random-float 100 < 70; 5 working days per week, so around 70% e 30%
                          ; go to work (if not already there, due to rule ($) above)
                          [if activateNursingHomes and [pcolor] of patch-here != orange
                              [
                               let destination one-of patches with [pcolor = orange]
                               set new-xcor [pxcor] of destination  - 0.3 + random-float 0.6
                               set new-ycor [pycor] of destination  - 0.3 + random-float 0.6
                               set xcor-1 xcor set ycor-1 ycor
                               setxy new-xcor new-ycor
                              ]
                          ]
                          ; move around to any place, but if already in a white or cyan + 3 place make short random movements
                          [ifelse [pcolor] of patch-here = orange

                           [let destination one-of patches with [pcolor = white or pcolor = cyan + 3]
                            setxy [pxcor] of destination - 0.3 + random-float 0.6
                                  [pycor] of destination - 0.3 + random-float 0.6]

                           [set iw 0
                            while [new-xcor < (min-pxcor - 0.499) or new-xcor > (max-pxcor + 0.499) or
                               new-ycor < (min-pycor - 0.499) or new-ycor > (max-pycor + 0.499)
                               or [off] of patch new-xcor new-ycor = True or [pcolor] of patch new-xcor new-ycor = yellow
                              ]
                               [set iw iw + 1 if iw > 100000 and inspace and [pcolor] of patch new-xcor new-ycor = yellow [makeJump xcor ycor]
                                ; the control above is to avoid situation of perpetual while which, sometime, occurred here, but also to preserve
                                ; long repetition cycles effectively finishing and so  determining the random values sequence (in a back-compatibility effort)
                                if inSpace and [off] of patch-here = True [makeJump xcor ycor]
                                set new-xcor xcor + random-normal 0 sdRandomMovement / sdReductionF set new-ycor ycor + random-normal 0 sdRandomMovement / sdReductionF]

                           set xcor-1 xcor set ycor-1 ycor
                           setxy new-xcor new-ycor]
                           ]
                          ]
                        ]


;----------------------------------------
                    ; simple defaulf case
                    if not chosen
                        [if temporaryDestination [set xcor xcor-1 set ycor ycor-1]

                        ; moving around
                        set chosen True
                        set iw 0
                        while [new-xcor < (min-pxcor - 0.499) or new-xcor > (max-pxcor + 0.499) or
                               new-ycor < (min-pycor - 0.499) or new-ycor > (max-pycor + 0.499)
                               or [off] of patch new-xcor new-ycor = True or [pcolor] of patch new-xcor new-ycor = yellow
                              ]
                               [set iw iw + 1 if iw > 100000 and inspace and [pcolor] of patch new-xcor new-ycor = yellow [makeJump xcor ycor]
                                ; the control above is to avoid situation of perpetual while which, sometime, occurred here, but also to preserve
                                ; long repetition cycles effectively finishing and so  determining the random values sequence (in a back-compatibility effort)
                                if inSpace and [off] of patch-here = True [makeJump xcor ycor]
                                set new-xcor xcor + random-normal 0 sdRandomMovement / sdReductionF set new-ycor ycor + random-normal 0 sdRandomMovement / sdReductionF]

                        ifelse temporaryDestination [set temporaryDestination False]
                                                    ; to switch off temporary status if coming from visiting NHs or going to school
                                                    [set xcor-1 xcor set ycor-1 ycor]
                        setxy new-xcor new-ycor
                      ]

                    ]
                 ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to evolveInfectedPeople

   ask people with [color = red]
                          [if deadTick = ticks [
                                                set newDeadCount newDeadCount + 1
                                                set cumulativeDeadCount cumulativeDeadCount + 1
                                                if inNursingHome [set cumulativeDeadCountNH cumulativeDeadCountNH + 1]
                                                die
                                               ]
                          ]

  ask people with [color = red]
      [if finishingInfection + 1 = ticks
                         [set color turquoise - 1 set size size / 1.2
                          set infectedInHospital False
                          set infectedAtHome     False

                          set newGreenCount newGreenCount + 1
                          set newGreenCount* newGreenCount* + 1
                          set cumulativeGreenCount cumulativeGreenCount + 1
                          set cumulativeGreenCount* cumulativeGreenCount* + 1
                          if inNursingHome [set cumulativeGreenCountNH cumulativeGreenCountNH + 1
                                            set cumulativeGreenCount*NH cumulativeGreenCount*NH + 1]]
      ]

    ask people with [color = violet]
      [if finishingInfection + 1 = ticks
                         [set color green + 1 set size size / 1.2
                          set infectedInHospital False
                          set infectedAtHome     False

                          set newGreenCount newGreenCount + 1
                          set cumulativeGreenCount cumulativeGreenCount + 1
                          if inNursingHome [set cumulativeGreenCountNH cumulativeGreenCountNH + 1]]
     ]


end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to activateDelayedInfections

  set countD 0
  repeat populationSize  [
                        if is-turtle? person countD [ask person countD [if delayedStart = ticks ; is-turtle? to bypass dead agents
                           [
                            ;show word "found " who
                            set delayedStart 0  ;clean delayed condition
                            ifelse pcolor != pink and pcolor != orange ; not in hospital or nursing home
                               [
                                set startingInfection ticks ; no incubation here
                                set finishingInfection delayedFinish
                                set delayedFinish 0

                                ;output seguence prepaparion (set below or ";" already set or ";;" to be set in haveInfectedPeopleCommingOut:
                                set myInfectingAgent 0                   ;external unindified agent
                                set myInfectingAgentContagion# 0         ;as above
                                set myInfectingAgentWhereInfectedColor 0 ;as above
                                ;who
                                ;myFragilityRate
                                set myContagion# contagion# set contagion# contagion# + 1
                                set whereInfectedColor 0                 ;as above
                                set contagionTick ticks
                                ;startingInfection
                                ;finishingInfection
                                ;;symptomaticAsymptomaticStatus
                               ]
                               [
                                set delayedFinish 0 ; clean also if in hospital or NH
                               ]
                           ]]
                           ]
                            set countD countD + 1
                           ]


end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to haveInfectedPeopleComingOut
  ask people with [startingInfection = ticks]
               [if shape = "cylinder" or shape = "star" or shape = "triangle" ; regular by construction
                [
                 ifelse random-float 100 < asymptomaticRegularInfected% or activeVaccination
                 [set color violet set size size * 1.2 if displayStatus = "HideGrayPeople" [set hidden? False]
                  set newVioletCount newVioletCount + 1
                  set cumulativeVioletCount cumulativeVioletCount + 1
                  if inNursingHome [set cumulativeVioletCountNH cumulativeVioletCountNH + 1]
                  if student [set cumulativeVioletCountStudents cumulativeVioletCountStudents + 1]
                  if teacher [set cumulativeVioletCountTeachers cumulativeVioletCountTeachers + 1]]
                 [set color red set size size * 1.2 if displayStatus = "HideGrayPeople" [set hidden? False]
                  set newRedCount newRedCount + 1
                  set cumulativeRedCount cumulativeRedCount + 1
                  if inNursingHome [set cumulativeRedCountNH cumulativeRedCountNH + 1]
                  if student [set cumulativeRedCountStudents cumulativeRedCountStudents + 1]
                  if teacher [set cumulativeRedCountTeachers cumulativeRedCountTeachers + 1]
                  if random-float 100 < dead%rateRegularPeople [set deadTick round startingInfection
                                              + (finishingInfection - startingInfection) / 2]
                  ]
                  set symptomaticAsymptomaticStatus 1 if color = violet [set symptomaticAsymptomaticStatus 2]
                 ]

                if shape = "x"                       ; fragile
                 [
                  ifelse random-float 100 < asymptomaticFragileInfected% or activeVaccination
                  [set color violet set size size * 1.2 if displayStatus = "HideGrayPeople" [set hidden? False]
                  set newVioletCount newVioletCount + 1
                  set cumulativeVioletCount cumulativeVioletCount + 1
                  if inNursingHome [set cumulativeVioletCountNH cumulativeVioletCountNH + 1]]
                 [set color red set size size * 1.2 if displayStatus = "HideGrayPeople" [set hidden? False]
                  set newRedCount newRedCount + 1
                  set cumulativeRedCount cumulativeRedCount + 1
                  if inNursingHome [set cumulativeRedCountNH cumulativeRedCountNH + 1]
                  if random-float 100 < dead%rateFragilePeople [set deadTick round startingInfection
                                              + (finishingInfection - startingInfection) / 2]
                 ]
                  set symptomaticAsymptomaticStatus 1 if color = violet [set symptomaticAsymptomaticStatus 2]
               ]
               set csvList lput (list myInfectingAgent myInfectingAgentContagion# myInfectingAgentWhereInfectedColor
                                      who myFragilityRate myContagion# whereInfectedColor contagionTick
                                      startingInfection finishingInfection symptomaticAsymptomaticStatus myGroupForStatistics) csvList
              ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to diffuseInfectionExcludingFactoriesSchools
  ; not at school or in a factory

  ; using here a sort to keep consistency with the web version
  ask people with [color = red and pcolor != yellow] ; '!= yellow' means 'not at school'
     [set infectingAgent self
      foreach sort other people in-radius radiusOfInfection with [startingInfection = 0]
        [aPerson ->
          ask aPerson
               [if random-float 1 < probabilityOfGettingInfection * intrinsicSusceptibility
                [
                 set startingInfection  ticks + 1 + incubationPeriod
                 set finishingInfection startingInfection + minInfectionDuration
                        + random (maxInfectionDuration - minInfectionDuration)

                 ;output seguence prepaparion (set below or ";" already set of ";;" to be set in haveInfectedPeopleCommingOut:
                 set myInfectingAgent [who] of infectingAgent
                 set myInfectingAgentContagion# [myContagion#] of infectingAgent
                 set myInfectingAgentWhereInfectedColor [whereInfectedColor] of infectingAgent
                 ;who
                 ;myFragilityRate
                 set myContagion# contagion# set contagion# contagion# + 1
                 set whereInfectedColor pcolor
                 set contagionTick ticks
                 ;startingInfection
                 ;finishingInfection
                 ;;symptomaticAsymptomaticStatus
                  ]
               ]
        ]
     ]

  ; using here a sort to keep consistency with the web version
  ask people with [color = violet and pcolor != yellow]  ; '!= yellow' means 'not at school'  ;more dangerous, IF people getting closer
      [set infectingAgent self
       foreach sort other people in-radius radiusOfInfection with [startingInfection = 0]
         [aPerson ->
           ask aPerson
               [if random-float 1 < probabilityOfGettingInfection * (1 + D% / 100) * intrinsicSusceptibility
                [
                 set startingInfection  ticks + 1 + incubationPeriod
                 set finishingInfection startingInfection + minInfectionDuration
                        + random (maxInfectionDuration - minInfectionDuration)

                 ;output seguence prepaparion (set below or ";" already set of ";;" to be set in haveInfectedPeopleCommingOut:
                 set myInfectingAgent [who] of infectingAgent
                 set myInfectingAgentContagion# [myContagion#] of infectingAgent
                 set myInfectingAgentWhereInfectedColor [whereInfectedColor] of infectingAgent
                 ;who
                 ;myFragilityRate
                 set myContagion# contagion# set contagion# contagion# + 1
                 set whereInfectedColor pcolor
                 set contagionTick ticks
                 ;startingInfection
                 ;finishingInfection
                 ;;symptomaticAsymptomaticStatus
                  ]
               ]
        ]
      ]
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to diffuseInfectionWithinFactories

  if cities&towns and activateFactories [
    ; red people never at work
    ask people with [color = violet and (pcolor = brown - 2 or pcolor = brown + 2)
                                    and openFactory]  ; brown means in factory, - 2 large, + 2 small
      [set infectingAgent self
       foreach sort other people-here with [startingInfection = 0]
        [aPerson ->
          ask aPerson
               [if random-float 1 < probabilityOfGettingInfection * (1 + D% / 100) * intrinsicSusceptibility ;more dangerous, people getting closer
                [set infectedAtWork True
                 set startingInfection  ticks + 1 + incubationPeriod
                 set finishingInfection startingInfection + minInfectionDuration
                        + random (maxInfectionDuration - minInfectionDuration)

                 ;output seguence prepaparion (set below or ";" already set of ";;" to be set in haveInfectedPeopleCommingOut:
                 set myInfectingAgent [who] of infectingAgent
                 set myInfectingAgentContagion# [myContagion#] of infectingAgent
                 set myInfectingAgentWhereInfectedColor [whereInfectedColor] of infectingAgent
                 ;who
                 ;myFragilityRate
                 set myContagion# contagion# set contagion# contagion# + 1
                 set whereInfectedColor pcolor
                 set contagionTick ticks
                 ;startingInfection
                 ;finishingInfection
                 ;;symptomaticAsymptomaticStatus
                  ]
               ]
        ]
      ]
  ]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to diffuseInfectionAtSchool

  if cities&towns and activateSchools [
    ; red people never at school
    ask people with [color = violet and pcolor = yellow]  ; yellow means in school
      [set infectingAgent self
       foreach sort other people-here with [startingInfection = 0 and classroom = [classroom] of myself]
        [aPerson ->
          ask aPerson
               [if random-float 1 < probabilityOfGettingInfection * (1 + D% / 100) * intrinsicSusceptibility ;more dangerous, people getting closer
                [set infectedAtSchool True
                 set startingInfection  ticks + 1 + incubationPeriod
                 set finishingInfection startingInfection + minInfectionDuration
                        + random (maxInfectionDuration - minInfectionDuration)

                 ;output seguence prepaparion (set below or ";" already set of ";;" to be set in haveInfectedPeopleCommingOut:
                 set myInfectingAgent [who] of infectingAgent
                 set myInfectingAgentContagion# [myContagion#] of infectingAgent
                 set myInfectingAgentWhereInfectedColor [whereInfectedColor] of infectingAgent
                 ;who
                 ;myFragilityRate
                 set myContagion# contagion# set contagion# contagion# + 1
                 set whereInfectedColor pcolor
                 set contagionTick ticks
                 ;startingInfection
                 ;finishingInfection
                 ;;symptomaticAsymptomaticStatus
                  ]
               ]
        ]
      ]
  ]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to diffuseInfectionAtHomeOrInNursingHomesOrInHospitalsNightStep
  ; in the same apartment on in the same nursing home or hospital, having a unique apartment

  ; people at home or in nursing home work in the same way automatically, people in hospitals are added,
  ; but look at the note [$$] below


  ; people infecting other people regularly in the apartment or in the same nursing home, in the nighty step of a day
  ; in the other steps other people in the houses can randomly meet
  ; about nursing homes, the infections to nursing home operators arise in the steps of the day in which they are
  ; there (the first one)

  ; using here a sort to keep consistency with the web version
  ask people with [color = red]
      [set infectingAgent self
       foreach sort other people with [xcorhome = [xcorhome] of myself and ycorhome = [ycorhome] of myself and
                              apartment = [apartment] of myself and startingInfection = 0]
        [aPerson ->
          ask aPerson
               [if random-float 1 < probabilityOfGettingInfection * intrinsicSusceptibility
                [
                 set startingInfection  ticks + 1 + incubationPeriod
                 set finishingInfection startingInfection + minInfectionDuration
                        + random (maxInfectionDuration - minInfectionDuration)

                 ;output seguence prepaparion (set below or ";" already set of ";;" to be set in haveInfectedPeopleCommingOut:
                 set myInfectingAgent [who] of infectingAgent
                 set myInfectingAgentContagion# [myContagion#] of infectingAgent
                 set myInfectingAgentWhereInfectedColor [whereInfectedColor] of infectingAgent
                 ;who
                 ;myFragilityRate
                 set myContagion# contagion# set contagion# contagion# + 1
                 set whereInfectedColor pcolor
                 set contagionTick ticks
                 ;startingInfection
                 ;finishingInfection
                 ;;symptomaticAsymptomaticStatus
                  ]
               ]
        ]
      ]

  ; using here a sort to keep consistency with the web version
  ask people with [color = violet]
      [set infectingAgent self
       foreach sort other people  with [xcorhome = [xcorhome] of myself and ycorhome = [ycorhome] of myself and
                               apartment = [apartment] of myself and startingInfection = 0]
        [aPerson ->
          ask aPerson
               [if random-float 1 < probabilityOfGettingInfection * (1 + D% / 100) * intrinsicSusceptibility ;more dangerous, people getting closer
                [
                 set startingInfection  ticks + 1 + incubationPeriod
                 set finishingInfection startingInfection + minInfectionDuration
                        + random (maxInfectionDuration - minInfectionDuration)

                 ;output seguence prepaparion (set below or ";" already set of ";;" to be set in haveInfectedPeopleCommingOut:
                 set myInfectingAgent [who] of infectingAgent
                 set myInfectingAgentContagion# [myContagion#] of infectingAgent
                 set myInfectingAgentWhereInfectedColor [whereInfectedColor] of infectingAgent
                 ;who
                 ;myFragilityRate
                 set myContagion# contagion# set contagion# contagion# + 1
                 set whereInfectedColor pcolor
                 set contagionTick ticks
                 ;startingInfection
                 ;finishingInfection
                 ;;symptomaticAsymptomaticStatus
                  ]
               ]
        ]
      ]

  ;case infected in hospital
  ;[$$] not relevant in this version of SIsaR as at night all resident person are red ones

  ; using here a sort to keep consistency with the web version
  ask people with [color = red and infectedInHospital]
      [ask other people-here with [startingInfection = 0]
               [if random-float 1 < probabilityOfGettingInfection * intrinsicSusceptibility
                [
                 set startingInfection  ticks + 1 + incubationPeriod
                 set finishingInfection startingInfection + minInfectionDuration
                        + random (maxInfectionDuration - minInfectionDuration) ]
               ]
      ]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to clearData
  set newRedCount   0
  set newVioletCount  0
  set newGreenCount 0
  set newGreenCount* 0
  set newDeadCount 0
  set cumulativeRedCount_t-1          cumulativeRedCount
  set cumulativeVioletCount_t-1       cumulativeVioletCount
  set cumulativeGreenCount_t-1        cumulativeGreenCount
  set cumulativeGreenCount_t-1*       cumulativeGreenCount*
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to makeCheckPoints

 if checkPoint1activeAt = ticks
  [set cp1cRed cumulativeRedCount
   set cp1cRedNH cumulativeRedCountNH
   set cp1cViolet cumulativeVioletCount
   set cp1cVioletNH cumulativeVioletCountNH
   set cp1cDead cumulativeDeadCount]

 if checkPoint2activeAt = ticks
  [set cp2cRed cumulativeRedCount
   set cp2cRedNH cumulativeRedCountNH
   set cp2cViolet cumulativeVioletCount
   set cp2cVioletNH cumulativeVioletCountNH
   set cp2cDead cumulativeDeadCount]

 if checkPoint3activeAt = ticks
  [set cp3cRed cumulativeRedCount
   set cp3cRedNH cumulativeRedCountNH
   set cp3cViolet cumulativeVioletCount
   set cp3cVioletNH cumulativeVioletCountNH
   set cp3cDead cumulativeDeadCount]

  if checkPoint4activeAt = ticks
  [set cp4cRed cumulativeRedCount
   set cp4cRedNH cumulativeRedCountNH
   set cp4cViolet cumulativeVioletCount
   set cp4cVioletNH cumulativeVioletCountNH
   set cp4cDead cumulativeDeadCount]

  if checkPoint5activeAt = ticks
  [set cp5cRed cumulativeRedCount
   set cp5cRedNH cumulativeRedCountNH
   set cp5cViolet cumulativeVioletCount
   set cp5cVioletNH cumulativeVioletCountNH
   set cp5cDead cumulativeDeadCount]

  if checkPoint6activeAt = ticks
  [set cp6cRed cumulativeRedCount
   set cp6cRedNH cumulativeRedCountNH
   set cp6cViolet cumulativeVioletCount
   set cp6cVioletNH cumulativeVioletCountNH
   set cp6cDead cumulativeDeadCount]
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to makePlots

  set-current-plot "stocks in each tick"
  set-current-plot-pen "InfectedSymptomatic"
  plot count turtles with [color = red]
  set-current-plot-pen "InfectedAsymptomatic"
  plot count turtles with [color = violet]
  set-current-plot-pen "Recovered"
  plot count turtles with [color = green + 1 or color = turquoise - 1]
  set-current-plot-pen "Recovered-ex-symptomatic"
  plot count turtles with [color = turquoise - 1]
  set-current-plot-pen "Susceptible"
  plot count turtles with [color = gray + 2]
  set-current-plot-pen "Susceptible-Fragile"
  plot count turtles with [color = gray + 2 and shape = "x"]
  set-current-plot-pen "Dead"
  plot cumulativeDeadCount


  set-current-plot "stocks in each tick without Susceptible"
  set-current-plot-pen "InfectedSymptomatic"
  plot count turtles with [color = red]
  set-current-plot-pen "InfectedAsymptomatic"
  plot count turtles with [color = violet]
  set-current-plot-pen "Recovered"
  plot count turtles with [color = green + 1 or color = turquoise - 1]
  set-current-plot-pen "Recovered-ex-symptomatic"
  plot count turtles with [color = turquoise - 1]
  set-current-plot-pen "Dead"
 plot cumulativeDeadCount


  set-current-plot "Infected stock"
  set-current-plot-pen "InfectedSymptomatic"
  plot count turtles with [color = red]
  set-current-plot-pen "InfectedAsymptomatic"
  plot count turtles with [color = violet]

  set-current-plot "new cases"
  set-current-plot-pen "InfectedSymptomatic"
  plot newRedCount
;  set-current-plot-pen "InfectedAsymptomatic"
;  plot newVioletCount
;  set-current-plot-pen "Recovered"
;  plot newGreenCount
;  set-current-plot-pen "Recovered-ex-symptomatic"
;  plot newGreenCount*


  set-current-plot "cumulative values"
  set-current-plot-pen "InfectedSymptomatic"
  plot cumulativeRedCount
  set-current-plot-pen "InfectedAsymptomatic"
  plot cumulativeVioletCount
  set-current-plot-pen "Recovered"
  plot cumulativeGreenCount
  set-current-plot-pen "Recovered-ex-symptomatic"
  plot cumulativeGreenCount*
  set-current-plot-pen "Dead"
  plot cumulativeDeadCount

;  set-current-plot "cumulative infection rate"
;  set-current-plot-pen "rate"
;  ifelse cumulativeRedCount_t-1 + cumulativeVioletCount_t-1 = 0
;                                    [plot-pen-up plot 0]
;                                    [plot-pen-down plot (newRedCount + newVioletCount) /
;                                                        (cumulativeRedCount_t-1 + cumulativeVioletCount_t-1)]

  set-current-plot "Rt at the end of each day"
  set-current-plot-pen "x-axis"
  plot 1
  set-current-plot-pen "Rt-1"
  if cumulativeRedCount_t-1 = 0 or  newGreenCount* = 0 [plot-pen-up plot 0]
  if cumulativeRedCount_t-1 > 0 and newGreenCount* > 0 [plot-pen-down plot newRedCount / newGreenCount*]
  set-current-plot-pen "Rt-2"
  if cumulativeRedCount_t-1 + cumulativeVioletCount_t-1 = 0 or  newGreenCount = 0 [plot-pen-up plot 0]
  if cumulativeRedCount_t-1 + cumulativeVioletCount_t-1 > 0 and newGreenCount > 0 [plot-pen-down plot (newRedCount
                                                                                  + newVioletCount) / newGreenCount]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to decodeScript

  set myScript script

  if myScript = "" [stop]

  set lll []

  ; check missing final \n
  if last myScript != " " and last myScript != "\n"
       [user-message "no final return (enter( key pressed" stop]

  while [is-number? position  "\n" myScript] ; \n is internally a single character
    [set myScript replace-item position  "\n" myScript  myScript " "]

  ; substitute multiple blanks with a unique one
  while [is-number? position  "  " myScript]
    [set myScript remove-item position  "  " myScript  myScript]

  ; how many rows?
  let myScript2 remove " " myScript
  set howManyRows (length myScript - length myScript2) / 3

  repeat howManyRows
  [
    set n 0
    while [ n < 3]
    [
      set s substring myScript 0 position " " myScript

      let ss s

      if n = 0 or n = 2 [set s read-from-string s]

      repeat 1 + length ss [set myScript but-first myScript]
      ; 1 + to eliminate right blank

      set lll lput s lll

      set n n + 1
    ]
  ]

  ; substrings
  set l []

  set nsub length lll / 3

  repeat nsub
  [ set sl []
    repeat 3 [

    set sl lput first lll sl   set lll but-first lll
             ]
    set l lput sl l
  ]

  set delayedStartList [] ; for the use in checkScript
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report fromSeedSequence [metaSeed]

  random-seed metaSeed
  let newSeed 0
  repeat run# [set newSeed random 100000000]
  report newSeed
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to checkScript

  set i 0 if ticks = 0 [set fileVacc 0]
  repeat nsub
    [
      if item 0 item i l = ticks [
                            if item 1 item i l = "seed"       [set seed item 2 item i l
                                                               if seed = -9999 [set seed int (timer * 1000)]
                                                               if seed != -9999 and seed < 0 [set seed fromSeedSequence abs seed]
                                                               if ticks = 0 [set seed0 seed]
                                                               ifelse ticks <= 1 [random-seed seed]
                                                                                 [user-message "use seed at  0 or 1" stop]]
                            if item 1 item i l = "seedD"      [set seedD item 2 item i l
                                                               if ticks != 0 [user-message "use seedD at  0" stop]]
                            if item 1 item i l = "vacc"       [set fileVacc item 2 item i l
                                                               if ticks != 0 [user-message "use vacc at  0" stop]]
                            if item 1 item i l = "add@"       [set delayedStartList lput (item 2 item i l) delayedStartList
                                                               if ticks != 0 [user-message "use add@ at  0" stop]]
                            if item 1 item i l = "csv"        [set choice item 2 item i l if choice = 0 [set saveCSVfile False]
                                                                                          if choice = 1 [set saveCSVfile True]
                                                                                          if choice != 0 and choice != 1 [user-message "use csv with 0 or 1 values" stop]]
                            if item 1 item i l = "myStop"     [set myStop item 2 item i l]
                            if item 1 item i l = "ratio"      [set ratioInitialInfected% item 2 item i l]
                            if item 1 item i l = "fP%"        [set fragilePeople% item 2 item i l]
                            if item 1 item i l = "yP%"        [set youngPeople% item 2 item i l]
                            if item 1 item i l = "d%R"        [set dead%rateRegularPeople item 2 item i l]
                            if item 1 item i l = "d%F"        [set dead%rateFragilePeople item 2 item i l]
                            if item 1 item i l = "minInf"     [set minInfectionDuration item 2 item i l]
                            if item 1 item i l = "maxInf"     [set maxInfectionDuration item 2 item i l]
                            if item 1 item i l = "incub"      [set incubationPeriod item 2 item i l]
                            if item 1 item i l = "#mov"       [set #movementCyclesPerTick item 2 item i l]
                            if item 1 item i l = "%PeopleAny" [set %PeopleAnyTypeNotSymptomaticLeavingHome item 2 item i l]
                            if item 1 item i l = "%PeopleNot" [set %PeopleNotFragileNotSymptomaticLeavingHome item 2 item i l]
                            if item 1 item i l = "%Fac"       [set %openFactoriesAlsoWhenLimitationsOn item 2 item i l]
                            if item 1 item i l = "%St"        [set %Students item 2 item i l]
                               if item 1 item i l = "limlock"    [set choice item 2 item i l
                                                                if choice = 0 [set activateLimitations/LockdownFromTick 0]
                                                                if choice = 1 [set activateLimitations/LockdownFromTick item 0 item i l]
                                                                if choice != 0 and choice != 1 [user-message "use lock with 0 or 1 values" stop]]
                            if item 1 item i l = "aBZ"        [set choice item 2 item i l
                                                                if choice = 0 [set activateBufferZones False]
                                                                if choice = 1 [set activateBufferZones True]
                                                                if choice != 0 and choice != 1 [user-message "use aBZ with 0 or 1 values" stop]]
                            if item 1 item i l = "radius"     [set radiusOfInfection item 2 item i l]
                            if item 1 item i l = "sd"         [set sdRandomMovement item 2 item i l]
                            if item 1 item i l = "prob"       [set probabilityOfGettingInfection item 2 item i l]
                            if item 1 item i l = "aH"         [set choice item 2 item i l if choice = 0 [set activateHospitals False]
                                                                                          if choice = 1 [set activateHospitals True]
                                                                                          if choice != 0 and choice != 1 [user-message "use aH with 0 or 1 values" stop]]
                            if item 1 item i l = "aNH"        [set choice item 2 item i l if choice = 0 [set activateNursingHomes False]
                                                                                          if choice = 1 [set activateNursingHomes True]
                                                                                          if choice != 0 and choice != 1 [user-message "use aNH with 0 or 1 values" stop]]
                            if item 1 item i l = "aHOp"       [set choice item 2 item i l if choice = 0 [set activateHOperators False]
                                                                                          if choice = 1 [set activateHOperators True]
                                                                                          if choice != 0 and choice != 1 [user-message "use aHOp with 0 or 1 values" stop]]
                            if item 1 item i l = "aNHOp"      [set choice item 2 item i l if choice = 0 [set activateNHOperators False]
                                                                                          if choice = 1 [set activateNHOperators True]
                                                                                          if choice != 0 and choice != 1 [user-message "use aNHOp with 0 or 1 values" stop]]
                            if item 1 item i l = "aSch"       [set choice item 2 item i l if choice = 0 [set activateSchools False]
                                                                                          if choice = 1 [set activateSchools True]
                                                                                          if choice != 0 and choice != 1 [user-message "use aSch with 0 or 1 values" stop]]
                            if item 1 item i l = "sFW"        [set choice item 2 item i l if choice = 0 [set stopFragileWorkers False]
                                                                                          if choice = 1 [set stopFragileWorkers True]
                                                                                          if choice != 0 and choice != 1 [user-message "use sFW with 0 or 1 values" stop]]
                            if item 1 item i l = "fWH"         [ifelse item 0 item i l = 0
                                                                [set choice item 2 item i l if choice = 0 [set fragileWorkersAtHome False]
                                                                                          if choice = 1 [set fragileWorkersAtHome True]
                                                                                          if choice != 0 and choice != 1 [user-message "use fWH with 0 or 1 values" stop]]
                                                                [user-message "use fWH 1 uniquely at tick 0, as a counterfactual action (look at Info)" stop]
                                                              ]
                            if item 1 item i l = "hMRs"       [ifelse item 0 item i l = 0
                                                               [
                                                                set howManyClassroomsInASchool item 2 item i l
                                                                if howManyClassroomsInASchool != 1 and howManyClassroomsInASchool != 2
                                                                    [user-message "use hMRs with 1 or 2 values" stop]
                                                               ]
                                                               [user-message "use fWH 1 uniquely at tick 0" stop]
                                                              ]
                            if item 1 item i l = "aFac"       [ifelse item 0 item i l = 0
                                                               [
                                                                set choice item 2 item i l if choice = 0 [set activateFactories False]
                                                                                          if choice = 1 [set activateFactories True]
                                                                                          if choice != 0 and choice != 1 [user-message "use aFac with 0 or 1 values" stop]
                                                               ]
                                                               [user-message "use aFac uniquely at tick 0" stop]
                                                              ]
                            if item 1 item i l = "PvNH"       [set choice item 2 item i l if choice = 0 [set peopleVisitingNHs False]
                                                                                          if choice = 1 [set peopleVisitingNHs True]
                                                                                          if choice != 0 and choice != 1 [user-message "use PvNH with 0 or 1 values" stop]]
                            if item 1 item i l = "assH"       [set assignmentRateInfectedToHospitals% item 2 item i l]


                            if item 1 item i l = "flash"      [set choice item 2 item i l
                                                               if choice = 1 [output-show (word "at " ticks " t.inf.sym.   " cumulativeRedCount)]
                                                               if choice = 2 [output-show (word "at " ticks " t.inf.asy.   " cumulativeVioletCount)]
                                                               if choice = 3 [output-show (word "at " ticks " t.inf.sym.NH " cumulativeRedCountNH)]
                                                               if choice = 4 [output-show (word "at " ticks " t.inf.asy.NH " cumulativeVioletCountNH)]
                                                               if choice = 5 [output-show (word "at " ticks " t.dec.       " cumulativeDeadCount)]
                                                               if choice != 1 and choice != 2 and choice != 3 and choice != 4 and choice != 5
                                                                                          [user-message "use flash with 1, 2, 3 values" stop]
                                                              ]

                            if item 1 item i l = "aCP"        [set choice item 2 item i l
                                                               if choice = 1 [set checkPoint1activeAt item 0 item i l]
                                                               if choice = 2 [set checkPoint2activeAt item 0 item i l]
                                                               if choice = 3 [set checkPoint3activeAt item 0 item i l]
                                                               if choice = 4 [set checkPoint4activeAt item 0 item i l]
                                                               if choice = 5 [set checkPoint5activeAt item 0 item i l]
                                                               if choice = 6 [set checkPoint6activeAt item 0 item i l]
                                                               if choice != 1 and choice != 2 and choice != 3 and choice != 4 and choice != 5 and choice != 6
                                                                             [user-message "use aCP with 1, 2, 3, 4 values" stop]
                                                              ]



                                 ]
      set i i + 1
    ]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to assignPeopleHome

  ; home in T or C / without T or C
  ifelse cities&towns
    [ifelse random-float 1 < 0.35 [set myBuilding one-of patches with [pcolor = cyan + 3 and plabel = "C"]]
                                  [set myBuilding one-of patches with [pcolor = cyan + 3 and plabel = "T"]]]
    [set myBuilding one-of patches with [pcolor = cyan + 3]]

    ; let specify the apartment in a buiding (patch)
    set apartment random 6 ; apartments with number 0-5

    set xcorHome [pxcor] of myBuilding
    set ycorHome [pycor] of myBuilding

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to movetoWhiteCyan
ifelse cities&towns

  [setxy xcorHome ycorHome
   set pa patch-here ; to esclude patch-here, avoiding (i) a wrong set of xcor-1 ycor-1 in makeMovement: (ii) a sintax error arising
                     ; if we use directly patch-here below
   set rr one-of sort radiusVector
   set pos one-of sort patches with [not (pxcor = [pxcor] of pa and
                                          pycor = [pycor] of pa)  and (pcolor = white or pcolor = cyan + 3)] in-radius rr
    move-to pos
    setxy xcor - 0.5 + random-float 0.99
          ycor - 0.5 + random-float 0.99
  ]
  ; sort above is ncessary to mantain consinstency betwewn the desk and the web versions of the program


  [move-to one-of patches with [pcolor = white or pcolor = cyan + 3]
                setxy xcor - 0.5 + random-float 0.99
                      ycor - 0.5 + random-float 0.99
  ]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to checkBufferZones

   set offPatches count patches with [off]

   if activateBufferZones
   [if offPatches = 0
    ;creates "off" zones
    [
     ask patches with [pcolor = orange or pcolor = pink][set off True]
     ask patches with [pcolor = orange or pcolor = pink][ask neighbors[
                                          if pcolor != cyan + 3
                                                    [set off True]]]]
  ]

  if not activateBufferZones
   [if offPatches > 0
    ;delete "off" zones
    [
        ask patches [set off False]
    ]
  ]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to makeJump [nx ny]
  ask patch nx ny [set openPatch one-of neighbors with [not off and pcolor != yellow]
                   if openPatch = nobody [set openPatch one-of neighbors with [pcolor != yellow]]
                  ]
  set xcor [pxcor] of openPatch
  set ycor [pycor] of openPatch

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report inSpace

  ifelse new-xcor >= (min-pxcor - 0.499) and new-xcor <= (max-pxcor + 0.499) and
         new-ycor >= (min-pycor - 0.499) and new-ycor <= (max-pycor + 0.499)
  [report True]
  [report False]

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setDisplay
  if displayStatus != displayChoice
  [set displayStatus displayChoice

   if displayChoice = "HidePeople"     [let j 0 repeat populationSize[if person j != nobody [ask person j [set hidden? True]] set j j + 1]]
   if displayChoice = "HideGrayPeople" [let j 0 repeat populationSize[if person j != nobody [ask person j [if color = gray + 2 [set hidden? True]]] set j j + 1]]
   if displayChoice = "DisplayPeople"  [let j 0 repeat populationSize[if person j != nobody [ask person j [set hidden? False]] set j j + 1]]
   if displayChoice = "SimpleMap"      [let k 0
                                        repeat howManyDecorativeAgents [ask aDecAg (populationSize + k) [set hidden?  True] set k k + 1]]
   if displayChoice = "DressingMap"    [let k 0
                                        repeat howManyDecorativeAgents [ask aDecAg (populationSize + k) [set hidden? False] set k k + 1]]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to makeFinalOutputAndCalculations

      output-show (word "finishing at " ticks)
      output-show (word "at " ticks " t.inf.sym.   " cumulativeRedCount)
      output-show (word "at " ticks " t.inf.asy.   " cumulativeVioletCount)
      output-show (word "at " ticks " t.inf.sym.NH " cumulativeRedCountNH)
      output-show (word "at " ticks " t.inf.asy.NH " cumulativeVioletCountNH)
      output-show (word "at " ticks " t.dec.       " cumulativeDeadCount)
      set TinfSch count people with [teacher and infectedAtSchool]
      set SinfSch count people with [student and infectedAtSchool]
      set TinfNotSch count people with [teacher and startingInfection != 0 and not infectedAtSchool]
      set SinfNotSch count people with [student and startingInfection != 0 and not infectedAtSchool]

      if saveCSVfile [set fileName user-input "enter the name of the csv file for sequential analysis, without extension" csv:to-file word fileName ".csv" csvList] ;comment for web version
      if saveCSVfile and fileVacc != 0 [csv:to-file word fileName "_list.csv" outList0] ;comment for web version

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to checkPeopleSize

  if peopleSize != peopleSize0
      [
       if peopleSize = "Big"
        [set z 0 repeat populationSize [set tmpAg person z if tmpAg != nobody [ask tmpAg [set size size * 2.2]] set z z + 1]]

       if peopleSize = "Regular"
        [set z 0 repeat populationSize [set tmpAg person z if tmpAg != nobody [ask tmpAg [set size size / 2.2]] set z z + 1]]

       set peopleSize0 peopleSize
      ]
end
@#$#@#$#@
GRAPHICS-WINDOW
188
162
655
630
-1
-1
9.0
1
10
1
1
1
0
0
0
1
0
50
0
50
1
1
1
ticks
30.0

SLIDER
13
55
185
88
populationSize
populationSize
0
10000
4350.0
10
1
NIL
HORIZONTAL

BUTTON
0
126
64
160
NIL
setup\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
192
107
386
160
colors and their meaning\npink, hospitals; orange, nursing homes; cyan, houses; yellow, schools; brown, fact.s-off.s-shops (small, large)
9
0.0
1

BUTTON
1
222
64
255
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
15
397
177
430
ratioInitialInfected%
ratioInitialInfected%
0
1
0.05
0.05
1
NIL
HORIZONTAL

TEXTBOX
798
73
996
168
gray - Susceptible (total)\ngray darker - Susceptible-Fragile\nred  - Infected symptomatic\nviolet  - Infected asymptomatic\nlime - Recovered (total)\nturquoise - Recovered-ex-symptomatic\nblack  - Dead
9
0.0
1

TEXTBOX
17
432
167
450
infected individuals infos
11
0.0
1

SLIDER
15
450
180
483
minInfectionDuration
minInfectionDuration
1
20
10.0
1
1
NIL
HORIZONTAL

SLIDER
15
488
179
521
maxInfectionDuration
maxInfectionDuration
minInfectionDuration + 1
35
25.0
1
1
NIL
HORIZONTAL

PLOT
657
163
964
430
stocks in each tick
t
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"InfectedSymptomatic" 1.0 0 -2674135 true "" ""
"Recovered" 1.0 0 -13840069 true "" ""
"Susceptible" 1.0 0 -4539718 true "" ""
"Susceptible-Fragile" 1.0 0 -9276814 true "" ""
"InfectedAsymptomatic" 1.0 0 -8630108 true "" ""
"Dead" 1.0 0 -16777216 true "" ""
"Recovered-ex-symptomatic" 1.0 0 -8862290 true "" ""

SLIDER
17
582
183
615
radiusOfInfection
radiusOfInfection
0
1
0.2
0.05
1
NIL
HORIZONTAL

PLOT
968
163
1236
431
new cases
t
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"InfectedSymptomatic" 1.0 0 -2674135 true "" ""
"Recovered" 1.0 0 -13840069 true "" ""
"InfectedAsymptomatic" 1.0 0 -8630108 true "" ""
"Recovered-ex-symptomatic" 1.0 0 -8862290 true "" ""

TEXTBOX
22
563
172
581
infection data\n
11
0.0
1

SLIDER
16
618
184
651
sdRandomMovement
sdRandomMovement
0
1
0.3
0.05
1
NIL
HORIZONTAL

PLOT
969
435
1238
633
cumulative values
t
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"InfectedSymptomatic" 1.0 0 -2674135 true "" ""
"Recovered" 1.0 0 -13840069 true "" ""
"InfectedAsymptomatic" 1.0 0 -8630108 true "" ""
"Dead" 1.0 0 -16777216 true "" ""
"Recovered-ex-symptomatic" 1.0 0 -8862290 true "" ""

SLIDER
16
525
180
558
incubationPeriod
incubationPeriod
0
15
7.0
1
1
NIL
HORIZONTAL

INPUTBOX
195
10
320
70
healthcareOperators%
1.0
1
0
Number

INPUTBOX
320
10
408
70
fragilePeople%
35.0
1
0
Number

TEXTBOX
1
10
214
66
S.I.s.a.R. v.0.9.6.1 - 20210731\nAuthors, COPYRIGHT and \nLICENSE in Info sheet
11
0.0
1

PLOT
659
435
965
631
stocks in each tick without Susceptible
NIL
NIL
0.0
10.0
0.0
0.3
true
false
"" ""
PENS
"Dead" 1.0 0 -16777216 true "" ""
"InfectedSymptomatic" 1.0 0 -2674135 true "" ""
"Recovered" 1.0 0 -13840069 true "" ""
"InfectedAsymptomatic" 1.0 0 -8630108 true "" ""
"Recovered-ex-symptomatic" 1.0 0 -8862290 true "" ""

INPUTBOX
713
10
877
70
probabilityOfGettingInfection
0.02
1
0
Number

PLOT
1657
13
1927
212
Rt at the end of each day
NIL
NIL
0.0
10.0
0.0
3.0
true
false
"" ""
PENS
"Rt-1" 1.0 0 -2674135 true "" ""
"x-axis" 1.0 0 -16777216 true "" ""
"Rt-2" 1.0 0 -6459832 true "" ""

TEXTBOX
1693
39
1888
84
basic reproduction ratios\n(horizontal black line indicates 1)\nlook at the info sheet for calculations
9
0.0
1

TEXTBOX
995
80
1171
162
The infection can generate a (i) symptomatic (red) or (ii) an asymptomatic (violet) subject, following the condition of (i) fragile or (ii) regular people
11
0.0
1

SWITCH
13
89
183
122
cities&towns
cities&towns
0
1
-1000

TEXTBOX
887
16
1083
68
asymptomatic Infected people modify the probability of getting the infection; let suppose an increase or a decrease of D%
10
0.0
1

INPUTBOX
1084
10
1134
70
D%
-50.0
1
0
Number

INPUTBOX
538
10
710
70
asymptomaticRegularInfected%
95.0
1
0
Number

INPUTBOX
538
72
707
132
asymptomaticFragileInfected%
20.0
1
0
Number

TEXTBOX
3
164
187
221
If cities&towns is On, after hitting 'setup', we have to wait a few seconds while the program sets the initial people distribution
10
115.0
1

MONITOR
1239
205
1373
250
InfectedAsymptomatic
count people with [color = violet]
17
1
11

MONITOR
1239
160
1374
205
InfectedSymptomatic
count people with [color = red]
17
1
11

MONITOR
1239
252
1374
297
Recovered
count people with [color = green]
17
1
11

PLOT
1239
300
1401
435
Infected stock
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"InfectedAsymptomatic" 1.0 1 -8630108 true "" ""
"InfectedSymptomatic" 1.0 1 -2674135 true "" ""

TEXTBOX
1269
325
1318
354
in each\ntick
11
0.0
1

TEXTBOX
127
209
183
340
Stops if no more \ninfected exist\nor /pause/\nat tick\nmyStop (if not 0)\n
11
0.0
1

TEXTBOX
394
71
529
128
Zoom \"Smaller\" to see more items\n\n
13
15.0
1

MONITOR
1250
438
1387
483
%Recov/TotPop
100 * cumulativeGreenCount / populationSize
17
1
11

TEXTBOX
1246
482
1434
505
also counting asymptomatic\n
11
0.0
1

MONITOR
1250
502
1387
547
%Recov/TotPop
100 * cumulativeGreenCount* / populationSize
17
1
11

TEXTBOX
1246
548
1434
571
not counting asymptomatic
11
0.0
1

SLIDER
16
654
184
687
dead%rateRegularPeople
dead%rateRegularPeople
0
10
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
14
689
185
722
dead%rateFragilePeople
dead%rateFragilePeople
0
10
3.0
0.1
1
NIL
HORIZONTAL

MONITOR
1250
566
1388
611
Deceased
cumulativeDeadCount
17
1
11

TEXTBOX
1245
613
1433
636
from symptomatic people\n
11
0.0
1

TEXTBOX
1159
121
1259
159
stocks at current tick
13
15.0
1

INPUTBOX
68
222
123
283
myStop
0.0
1
0
Number

SLIDER
197
72
388
105
#movementCyclesPerTick
#movementCyclesPerTick
1
10
2.0
1
1
NIL
HORIZONTAL

SLIDER
623
638
897
671
activateLimitations/LockdownFromTick
activateLimitations/LockdownFromTick
0
100
20.0
1
1
NIL
HORIZONTAL

TEXTBOX
657
676
1042
696
0 means NEVER; > 0, from that tick
13
15.0
1

BUTTON
665
773
887
808
NIL
hardFinishLockdown
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
1022
697
1271
757
%PeopleAnyTypeNotSymptomaticLeavingHome
80.0
1
0
Number

TEXTBOX
955
678
1020
746
while the lockdown is on
13
15.0
1

INPUTBOX
1022
636
1270
696
%PeopleNotFragileNotSymptomaticLeavingHome
90.0
1
0
Number

TEXTBOX
722
79
792
142
Scroll on the right to prepare a script
12
15.0
1

TEXTBOX
1412
13
1602
293
HOW TO PREPARE A SCRIPT\nInsert the sequences of three items in the window on the right, by row as in the example here. You also can use the script to set the initial values (and keep note of them). Do not forget to end the last row with the return (o enter) key.\n\ntick shortName value\n60 %PeopleAny 20\n70 %PeopleNot 40\n...   \n(you can use empty rows to improve readability)
12
0.0
1

TEXTBOX
1674
227
1910
295
HOW TO SAVE THE INTERFACE\t\nRight click in an emply space and select Export interface, then choose file name and folder
12
115.0
1

TEXTBOX
1674
308
2049
1105
LONG NAMES / SHORT NAMES / VALUES\n\nvaccination vacc name (without quotation signs)\nrandom-seed seed anIntegerNumber \nyou can set the seed only at tick 0 or at tick 1\n[special setting: value -9999 or other negative v.\nto have the system choosing internally the seed]\nsaveCSVfile csv  (0 off, 1 on)\nmyStop    myStop 1...1000\n\ndead%rateRegularPeople d%R 0...10\ndead%rateFragilePeople d%F 0...10\nfragilePeople% fP% 0...100\nyoungPeople% yP% 0...100\nratioInitialInfected% ratio 0...100\nminInfectionDuration minInf 1...20\nmaxInfectionDuration maxInf minInf...35\nincubationPeriod incub 0...15\n#movementCyclesPerTick #mov 1...10\n\n%PeopleAnyTypeNotSymptomaticLeavingHome %PeopleAny 1...100\n%PeopleNotFragileNotSymptomaticLeavingHome %PeopleNot 1...100\n\nactivateHospitals aH  (0 off, 1 on)\nactivateNursingHomes aNH  (0 off, 1 on)\nactivateHOperators aHOp  (0 off, 1 on)\nactivateNHOperators aNHOp  (0 off, 1 on)\nactivateSchools aSch  (0 off, 1 on)\n%Students %St (0..100)\nhowManyRoomsInEachSchool hMRs (1 or 2)\nactivateFactories aFac  (0 off, 1 on)\n%openFactoriesWhenLimitationsOn %Fac 0...100\nfragileWorkersAtHome fWH (0 off, 1 on)\nstopFragileWorkers sFW (0 off, 1 on)\n\nlimitations/Lockdown limlock (0 off, 1 on)\nactivateBufferZones aBZ (0 off, 1 on)\n\nradiusOfInfection radius 0...1\nsdRandomMovement sd 0...1\nprobabilityOfGettingInfection prob 0...1\n\npeopleVisitingNHs pvNH (0 off, 1 on)\nassignmentRateInfectedToHospitals% assH 0..100\n\nflash flash (1 t.inf.sym., 2 t.inf.asy., 3  t.inf.sym.NH, 4  t.inf.asy.NH, 5 t.dec.)\nactivateCheckPoint aCP  (1 first, 2 second, 3 third, 4 fourth)\n
12
0.0
1

SWITCH
219
634
381
667
activateHospitals
activateHospitals
0
1
-1000

SWITCH
427
672
610
705
activateNursingHomes
activateNursingHomes
0
1
-1000

INPUTBOX
1403
305
1663
1108
script
0 seed  66507108\n1 seed   94261776\n0 seedD  44424105\n0 add@ 211\n0 add@ 211\n0 csv 0\n0 d%R 1\n0 d%F 3\n0 limlock 0\n0 myStop 0\n0 ratio 0.05\n0 minInf 10\n0 maxInf 25\n0 incub 7\n0 #mov 3\n0 radius 0.20\n0 sd 0.3\n0 prob 0.05\n\n0 aFac 1\n0 %Fac 100\n0 %PeopleNot 0\n0 %PeopleAny 0\n0 aSch 1\n0 %St 100\n0 hMRs 1\n0 fP% 35\n0 yP% 21\n0 aBZ 0\n0 aH 1\n0 aNH 1\n0 aHOp 1\n0 aNHOp 1\n0 fWH 0\n0 PvNH 1\n0 assH 100\n133 assH 33\n0 sFW 0\n\n1 aSch 1\n17 aSch 0\n\n20 limlock 1\n20 %PeopleAny 90\n28 %PeopleAny 80\n\n31 %PeopleAny 0\n31 %PeopleNot 80\n31 aBZ 1\n\n35 %PeopleNot 70\n35 sd 0.15\n35 #mov 2\n\n36 %PeopleNot 65\n\n38 PvNH 0\n38 %PeopleNot 15\n38 #mov 1\n38 %Fac 40\n\n42 %PeopleNot 25\n\n49 prob 0.02\n49 %Fac 20\n\n84 %PeopleNot 30\n84 %Fac 70\n\n106 %PeopleNot 0\n106 %PeopleAny 80\n106 sd 0.30\n106 #mov 2\n106 %Fac 100\n\n110 %PeopleAny 95\n110 #mov 3\n112 %PeopleAny 85\n112 #mov 2\n\n117 %PeopleAny 95\n117 #mov 3\n121 %PeopleAny 90\n121 #mov 2\n122 #mov 3\n149 prob 0.035\n\n104 aSch 0\n\n131 aSch 0\n\n225 aSch 1\n\n253 #mov 2\n\n259 %PeopleAny 90\n\n266 %PeopleAny 80\n266 %Fac 90\n266 prob 0.02\n\n277 %PeopleAny 50\n277 %Fac 70\n277 #mov 1\n277 sd 0.15\n277 PvNH 0\n\n302 %Fac 80\n302 %PeopleAny 70\n302 %PeopleNot 90\n302 #mov 2\n\n277 %St 50\n\n314 #mov 3\n\n325 aSch 0\n\n320 %PeopleAny 90\n320 %Fac 90\n320 #mov 3\n\n322 %Fac 90\n322 #mov 2\n\n325 %Fac 30\n325 %PeopleAny 50\n325 %PeopleNot 50\n325 #mov 1\n\n329 %Fac 90\n329 #mov 2\n329 %PeopleAny 80\n\n332 %Fac 30\n332 %PeopleAny 50\n332 %PeopleNot 50\n332 #mov 1\n\n336 %Fac 90\n336 #mov 2\n336 %PeopleAny 80\n\n337 %Fac 30\n337 %PeopleAny 50\n337 %PeopleNot 50\n337 #mov 1\n\n339 PvNH 1\n339 %PeopleAny 80\n339 %PeopleNot 100\n339 #mov 3\n339 %Fac 100\n339 aSch 1\n339 %St 50\n\n349 %PeopleNot 90\n349 #mov 2\n349 sd 0.30\n\n350 %St 50\n\n\n38 flash 1\n38 flash 2\n38 flash 3\n38 flash 4\n38 flash 5\n106 flash 1\n106 flash 2\n106 flash 3\n106 flash 4\n106 flash 5\n137 flash 1\n137 flash 2\n137 flash 3\n137 flash 4\n137 flash 5\n167 flash 1\n167 flash 2\n167 flash 3\n167 flash 4\n167 flash 5\n198 flash 1\n198 flash 2\n198 flash 3\n198 flash 4\n198 flash 5\n\n229 flash 1\n229 flash 2\n229 flash 3\n229 flash 4\n229 flash 5\n\n413 flash 1\n413 flash 2\n413 flash 3\n413 flash 4\n413 flash 5\n\n30   aCP 1\n119 aCP 2\n230 aCP 3\n316 aCP 4\n364 aCP 5\n453 aCP 6\n
1
1
String

TEXTBOX
70
125
179
167
About the seed of the random numbers see Info
9
0.0
1

MONITOR
1190
70
1297
115
Infected in n. h.
count people with [color = red and pcolor = orange]
17
1
11

MONITOR
1239
114
1361
159
Infected fragile
count people with [color = red and shape = \"x\"]
17
1
11

SWITCH
427
745
612
778
peopleVisitingNHs
peopleVisitingNHs
0
1
-1000

TEXTBOX
1282
675
1386
736
>>> up on the right R0/Rt values
14
15.0
1

INPUTBOX
220
718
410
778
assignmentRateInfectedToHospitals%
33.0
1
0
Number

TEXTBOX
1408
284
1646
303
>>> dictionary on the right
15
15.0
1

TEXTBOX
386
147
785
170
x, fragile persons; *, hospital operators; ∆, nursing home operators
11
0.0
1

SWITCH
220
672
380
705
activateHOperators
activateHOperators
0
1
-1000

SWITCH
427
708
609
741
activateNHoperators
activateNHoperators
0
1
-1000

MONITOR
1298
70
1392
115
infected in h,
count people with [color = red and pcolor = pink]
17
1
11

SWITCH
667
698
887
731
activateBufferZones
activateBufferZones
0
1
-1000

INPUTBOX
1159
10
1336
70
intrinsicSusceptibilityFactor
5.0
1
0
Number

TEXTBOX
1264
41
1333
59
look at Info
11
15.0
1

OUTPUT
2003
36
2333
412
12

TEXTBOX
1973
20
2064
38
Flash reports
11
0.0
1

INPUTBOX
409
10
514
70
youngPeople%
21.0
1
0
Number

SWITCH
13
729
196
762
activateSchools
activateSchools
0
1
-1000

CHOOSER
12
766
197
811
howManyClassroomsInASchool
howManyClassroomsInASchool
1 2
0

MONITOR
1958
439
2064
484
seed (last one)
seed
0
1
11

TEXTBOX
2103
417
2253
445
Cumulative values \nfrom nursing homes
11
0.0
1

MONITOR
2101
642
2268
687
Deceased
cumulativeDeadCountNH
17
1
11

MONITOR
2101
446
2262
491
Infected symptomatic
cumulativeRedCountNH
17
1
11

MONITOR
2100
494
2264
539
Infected asymptomatic
cumulativeVioletCountNH
17
1
11

MONITOR
2100
542
2265
587
Recovered ex-symptomatoc
cumulativeGreenCount*NH
17
1
11

MONITOR
2101
591
2267
636
Recovered ex-asymptomatic
cumulativeGreenCountNH - cumulativeGreenCount*NH
17
1
11

TEXTBOX
2100
690
2250
718
Cumulative values \nfrom teachers and students
11
0.0
1

MONITOR
2099
720
2245
765
Sympomatic students
cumulativeRedCountStudents
17
1
11

MONITOR
2100
770
2245
815
Asympomatic students
cumulativeVioletCountStudents
17
1
11

MONITOR
2100
820
2245
865
Sympomatic teachers
cumulativeRedCountTeachers
17
1
11

MONITOR
2100
870
2246
915
Aympomatic teachers
cumulativeVioletCountTeachers
17
1
11

MONITOR
2070
922
2174
967
teachs inf. sch.
count people with [teacher and infectedAtSchool]
17
1
11

MONITOR
2177
923
2276
968
studs inf. sch.
count people with [student and infectedAtSchool]
17
1
11

CHOOSER
3
341
121
386
displayChoice
displayChoice
"HidePeople" "HideGrayPeople" "DisplayPeople" "SimpleMap" "DressingMap"
2

SWITCH
214
855
409
888
activateFactories
activateFactories
0
1
-1000

INPUTBOX
1022
774
1271
834
%openFactoriesAlsoWhenLimitationsOn
100.0
1
0
Number

SWITCH
215
815
409
848
fragileWorkersAtHome
fragileWorkersAtHome
1
1
-1000

CHOOSER
4
293
96
338
peopleSize
peopleSize
"Big" "Regular"
1

TEXTBOX
385
637
593
671
Zoom \"Larger\" to see more agents
13
15.0
1

TEXTBOX
386
106
539
148
red - symp,; violet - asymp.; turquoise - symp. rec.; green - asym. rec.
11
0.0
1

SWITCH
427
814
612
847
stopFragileWorkers
stopFragileWorkers
1
1
-1000

SLIDER
10
908
183
941
run#
run#
1
10000
1.0
1
1
NIL
HORIZONTAL

TEXTBOX
10
875
225
904
The number identifying a run in BehaviorSpace tool
11
0.0
1

SLIDER
13
815
197
848
%Students
%Students
0
100
50.0
1
1
NIL
HORIZONTAL

SLIDER
2410
97
2583
130
g1t1
g1t1
0
1
0.1
0.01
1
NIL
HORIZONTAL

TEXTBOX
2417
37
2605
65
quotas: slides used by GA, or by ourselves, if \"y\" in the setting file
11
0.0
1

SLIDER
2409
139
2582
172
g2t1
g2t1
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2411
181
2583
214
g3t1
g3t1
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2411
223
2583
256
g4t1
g4t1
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2409
267
2581
300
g5t1
g5t1
0
1
0.43
0.01
1
NIL
HORIZONTAL

SLIDER
2409
311
2581
344
g6t1
g6t1
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2409
351
2581
384
g7t1
g7t1
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2598
98
2770
131
g1t2
g1t2
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2790
100
2962
133
g1t3
g1t3
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2976
100
3148
133
g1t4
g1t4
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
3162
100
3334
133
g1t5
g1t5
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2599
140
2772
173
g2t2
g2t2
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2788
140
2961
173
g2t3
g2t3
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2975
142
3148
175
g2t4
g2t4
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
3162
142
3334
175
g2t5
g2t5
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2599
181
2771
214
g3t2
g3t2
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2787
181
2959
214
g3t3
g3t3
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2975
183
3147
216
g3t4
g3t4
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
3162
183
3334
216
g3t5
g3t5
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2599
222
2771
255
g4t2
g4t2
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2784
224
2956
257
g4t3
g4t3
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2975
224
3147
257
g4t4
g4t4
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
3161
224
3333
257
g4t5
g4t5
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2600
268
2772
301
g5t2
g5t2
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2784
267
2956
300
g5t3
g5t3
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2971
267
3143
300
g5t4
g5t4
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
3158
267
3330
300
g5t5
g5t5
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2597
312
2769
345
g6t2
g6t2
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2784
313
2956
346
g6t3
g6t3
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2970
313
3142
346
g6t4
g6t4
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
3154
314
3326
347
g6t5
g6t5
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2595
353
2767
386
g7t2
g7t2
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2782
355
2954
388
g7t3
g7t3
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
2967
357
3139
390
g7t4
g7t4
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
3152
357
3324
390
g7t5
g7t5
0
1
0.1
0.01
1
NIL
HORIZONTAL

BUTTON
2410
422
2509
456
hard finish
set theEnd true
NIL
1
T
OBSERVER
NIL
H
NIL
NIL
1

MONITOR
2415
495
2777
540
iniltal susceptible persons
susceptibleList0
17
1
11

MONITOR
2414
559
2777
604
susceptible p. when vaccination starts
susceptibleListV
17
1
11

BUTTON
2829
496
3165
531
LaTeX show initial susceptible in commadn center
show (word \" Susc. at t = 0 & \"\nitem 1 susceptibleList0 \" & \"\nitem 3 susceptibleList0 \" & \"\nitem 5 susceptibleList0 \" & \"\nitem 7 susceptibleList0 \" & \"\nitem 9 susceptibleList0 \" & \"\nitem 11 susceptibleList0 \" & \"\nitem 13 susceptibleList0 \" \\\\ \"\n)
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
2828
560
3237
595
LaTeX show susceptible when vacc. starts in command center
if ticks >= vaccinationDay \n[show (word \" Susc. when\\\\vacc. starts & \" \nitem 1 susceptibleListV \" & \"\nitem 3 susceptibleListV \" & \"\nitem 5 susceptibleListV \" & \"\nitem 7 susceptibleListV \" & \"\nitem 9 susceptibleListV \" & \"\nitem 11 susceptibleListV \" & \"\nitem 13 susceptibleListV \" \\\\ \"\n)]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
2828
444
3164
477
db show initial susceptibe in command center
show (word \"[ \"\nitem 1 susceptibleList0 \", \"\nitem 3 susceptibleList0 \", \"\nitem 5 susceptibleList0 \", \"\nitem 7 susceptibleList0 \", \"\nitem 9 susceptibleList0 \", \"\nitem 11 susceptibleList0 \", \"\nitem 13 susceptibleList0 \" ]\"\n)
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## WHY THE NAME?

    S - Susceptible
    I - Infected
    s - symptomatic
    a - asymptomatic
    R - Recovered

with capital letters referring the *classic* [S.I.R. model](https://www.maa.org/press/periodicals/loci/joma/the-sir-model-for-spread-of-disease-the-differential-equation-model), also in [Wikipedia](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology).

#### Website

The model has a [website](https://terna.to.it/simul/SIsaR.html).

## PRELIMINARY NOTE

This is a simulation with random events, please do not take it as a sure forecasting machine: it is a reasoning machine, a sort of very complex "what if" mental experiment.

The New York Times offers us an analysis  on the [The Covid-19 Riddle: Why Does the Virus Wallop Some Places and Spare Others?](https://www.nytimes.com/2020/05/03/world/asia/coronavirus-spread-where-why.html) (open link). At the end of the article, we read: *Roll of the Dice - Finally, most experts agree that there may be no single reason for some countries to be hit and others missed. The answer is likely to be some combination of the above factors, as well as one other mentioned by researchers: sheer luck*.

In the same way, in the simulations run with this model we can have very different outcomes as we change the initial seed of the random numbers. Those values are determining mainly the movements at a tiny scale and so the interactions-infections chains.

We can use the model in a comparative way, observing different range of results with different initial conditions (parameters).

Finally, to have a reference at an actual situation, the model is related to the Piedmont scale, with 4,350 agents vs. 4.35 millions of inhabitants. The scale 1 to 1000 is over-represented in the case of schools, with their classrooms with a realistic number of students, apartments with a realistic quantity of inhabitants, and likewise workspaces, hospitals, nursing homes.

We look also to the time series of the total infected people in Piedmont.

![Total infections in Piedmont](https://terna.to.it/simul/andamento900.jpg)

## WEB VERSION

To generate the web version of the program, please comment the two row of code concluding with ";comment for web version".

## THINGS TO DO

Experiment with social distance (*radiusOfInfection*) considering this [review](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(20)31142-9/fulltext). A radius of 0.2 (the current one) is 20 meters, so we have to better calibrate with movement and probability, this is a key step in future 2.0 version. 

Better tune the death rate.

Verify the effects of different ways of behaving after the limitations/Lockdown.

Dilution effect and infection spread / there must be a definite viral load under which there is no infection, see on Nature, a paper wuth the virologist at Charité Universitätsmedizin Berlin, Christian Drosten [(Wölfel et al. 2020)] (https://www.nature.com/articles/s41586-020-2196-x).

Many viruses competing with different levels of danger: if we meet many infected persons, we have a higher probability of getting the worst one,

The COST of segregation, how much? and if we only segregate fragile people and extra fragile ones?

u - add unsusceptible people (which color?), are they necessary?

## VIEWING THE MODEL

In the desk version, we can use both the *continuous* view, observing all the agents' movements and the *on tick" one, faster but updating the screen only at the end of each tick. To set the view, we use a chooser in the top part of the screen.

In the web version we have uniquely the *on tick* feature.

The desk version also allows 3D view: right click on the map of the world and choose *Switch to 3D View*.

## HOW IT WORKS

### Scripting capability

The code can manage a script to set the parameter modifications occurring while the simulation is running. Explanations in the right side of the interface.

The simulation starts at tick 1, but we can already set the initial values at tick 0.

**A trick**: (i) via **setup** set all the default values for the experiments at tick 0 and then (ii) modify those that we want change before hitting **go**.

The scripting dictionary is reported in the Interface, close to the script window.

### Special item "add@"

With the "add@" item, mandatorily set at time 0, we add an infected agent (the choice is random, such as the symptomatic or asymptomatic type) at the day indicated after the item, e.g.,
0 add@ 211
means that at tick 211, one of the agents not in a hospital and not in a nursing home becomes infected; here is not relevant if the same agent was already infected in the past.

Using add@ it is mandatory to introduce in the script sequence 
0 seedD 'some value'

Here we have a third seed to avoid, until it is possible, the interferences with the random sequence of the setup and of the run. When the delayed infections emerge, the sequence of the run random numbers changes.
About 'some values': a value > 0 is a fixed seed; the value -9999 ask the system to generate randomly the seed; a negative non -9999 value asks the system to use that seed to generate a sequence of seedD in different runs of the simulation (sse: Special values: < 0 and ≠ -9999),

### Special item "flash"

With the "flash" item, followed by 1 or 2 or 3, we obtain a flash output in the window to the far right, with data at the beginning of the tick for: (1) total infected symptomatic people; (2) total infected asymptomatic people; (3) total infected symptomatic people in NH; (4) total infected asymptomatic people in NH;(5) total deceased.

### Special item "activate check point"

With the "aCP" (activate checkpoint) item, followed by n with n in [1,6] as [closed interval](https://abstractmath.org/MM/MMSetSpecific.htm#interval), we collect data at the tick of the aCP command. In BehaviorSpace, we can send (or not) the checkpoint data to the table of the results; we have six possible checkpoints, each reporting the cumulative numbers of infected symptomatic people, infected asymptomatic people, deceased people, related to the epidemics still alive at that tick. The data are at the beginning of the tick. 

We also have the measure of the final cumulative data that we can use in BehaviorSpace in the same way.

### The script in the model online

The script contained in the version of the model, that we find in the [repository](https://terna.to.it/simul/SIsaR.html), is related to a sequence of events and decisions referred to the  Piedmont region in Italy; day 1 there is February 3rd, 2020.

Day 0 is the setup time.

The epidemic simulated here could have a unique wave, concluding at day 105, but two new arrivals at days 211 restart the epidemic, concluding at 341; adding a third new arrival with '0 add@ 211', the epidemic would have a second huge wave, concluding at day 468. The lesson to learn is that of the enormous variability intrinsic to this kind of phenomenon.

### Random number seed

Using the scripting capability, we can choose a start for the sequence of the random numbers, i.e., the so called *seed*) that the program uses running the simulation.

Examples:

    0 seed 12345

or

    1 seed 1212113


If the referring tick is 0, we start using that sequence from the setup (both that seed and the successive sequence determine the world). If it is 1, we start using that sequence after the setup, and only the operating sequence is related to that seed. We can use both the definitions, one of them, none.

As ticks, we cannot use here values different from 0 or 1.

If we do not set the seed, NetLogo sets a value from a system's internal variable; the seed will be different in each execution. Unfortunately, we cannot access it to repeat the sequence if we need to do that.

#### Special value -9999

This value has a conventional meaning here: setting the seed of the random values to a casually internally generated seed, always different. In this case, we know the seed, and we can save it for future reuse.

At tick 0 the result is the same of not setting the seed (in case, NetLogo creates one internally).

An interesting use is that of setting a defined value at tick 0 and setting -9999 at tick 1. Example:

    0 seed 12345
    1 seed -9999

The setup, with the world and population creation, will be always the same; the run will be different. This choice can be useful for a comparative use of the model and to run a series of structurally identical experiments--with different random values--to obtain a distribution of results.

Useful also with the Tool *Behavior Space*.

#### Special values: < 0 and ≠ -9999

If we have

    0 seed v

or 

    1 seed v

with v < 0 and v ≠ -9999

the instrumental seed 'abs v' starts a sequence of random numbers in the interval 'random 100000000' where we extract the nth with n = run# (see below).


#### Using the Tool *BehaviorSpace*

We have two ways to obtain run replications in *BehaviorSpace*:

a) setting its internal variable Repetition to a value > 1 (with a flag in the related box 'Run combinations in sequential order');

b) setting its internal variable Repetition to 1 (the flag above is not relevant in this case) and placing in the first box ('Vary variables as follows) the reference to our instrumental variable *run#* (see below).

### Time management

The simulation stops when no more infected people (asymptomatic + symptomatic) exist.

If myStop is set to *a non-0 value*, the simulation stops at that number of ticks. Then we can set a new value and restart and the simulation will stop again. We can also pause the simulation hitting the *go* button when it is black and then hitting it again to restart.

### The flow of the model (*go* procedure)

The name of the procedures that we call for each *go* execution, so in a tick or *day* are the following (with self-explicative names).

First of all, with *checkBufferZones* we control whether to limit the movements of people in the zone of hospitals and nursing homes (same patch and conterminous patches), not considering those containing houses; the limitation is not operating for people directed to hospitals and nursing homes, as operators.

    checklimitations/Lockdown
    checkBufferZones

At the beginning of each tick, it is possible to modify the size of the agents and the display: look at "Hiding agents and map" below.

    checkPeopleSize
    setDisplay


The "day" starts at midnight, sending people at home, with *sendPeopleHome*, if they are around (in case, they will return in the same place where they were, in the "morning" of the same day).

Being people at home, or in hospital or in nursing home, *diffuseInfectionHomeOrInNursingHomesOrInHospitalsNightStep* produces infections with people within the same apartment, or within the same hospital or nursing home, for people who are resident there, always on a probabilistic basis.

    sendPeopleHome
    diffuseInfectionAtHomeOrInNursingHomesOrInHospitalsNightStep

we set to 0 the counter of the repetition of the actions (some of them are not repeated) within a day

    set countRepetitions 0

Then we start the repetitions, executed *#movementCyclesPerTick* times. 

Please note that the repetition in each tick applies to the steps *makeMovement* and *diffuseInfectionExcludingFactoriesSchools *, while we never repeat in the same tick the steps *evolveInfectedPeople* and *haveInfectedPeopleComingOut*, because they reflect the time development, so are linked to the tick value.

*makeMovement* simulates the movement of people, with or without lockout, with exceptions, following the *if* condition that we can examine in the beginning of this procedure. People at home do not move.

*evolveInfectedPeople* moves the subjects to the Recovered status, if the period as Infected is expired. In this procedure we control also the event fo the death.

*haveInfectedPeopleComingOut* makes explicit the subjects infected in previous periods and that have finished incubation.

*diffuseInfectionAtSchool* produces the spread of the infection for people within schools, if *activateSchool* in *on*; schools are only operating when countRepetitions is 0 (see below schools).

*diffuseInfectionWithinFActories* produces the spread of the infection for people within factories, if *activateSchool* in *on* and creates them; factories  are operating in the quota set by *%openFactoriesAlsoWhenLimitationsOn* and when countRepetitions is 0 (see below schools).


*diffuseInfectionExcludingFactoriesSchools* produces the spread of the infection (see below).

     while [countRepetitions < #movementCyclesPerTick]
     [

      makeMovement

     if countRepetitions = 0
      [

       evolveInfectedPeople
       haveInfectedPeopleComingOut
       diffuseInfectionAtSchool
       diffuseInfectionWithinFactories

      ]

      diffuseInfectionExcludingFactoriesSchools 

      set countRepetitions countRepetitions + 1
     ]

### Hiding agents and map

The chooser *setDisplay* hides or show the people or the map in the world, without effects on the simulation. If the run already finished, to activate the effects, we must hit *go*. 

Combining hiding of gray agent and the choice of size "Big" in the chooser *peopleSize* we can follow the evolution of the infection in the map.


### The calibration of the model

We calibrated the model using the data below, with specific references to Piedmont info. We shaped the behavior of the agents using a calendar of the events and accounting for limitation, lockdown, etc., also considering mobility data.

[Here](https://terna.to.it/simul/mobilityPiedmont.pdf) you have the mobility data (from Feb 15th, 2002, to May 29th, 2020), coming from Google LLC "Google COVID-19 Community Mobility Reports". https://www.google.com/covid19/mobility/ Accessed: Jun 6th, 2020.

### The infection diffusion (contagion)

With *ratioInitialInfected%* set to 0.05% the number of initial infected is a number close to 2, as 0.05% * 4350 => about 2. See below for the population size. Initially infected people bypass the incubation period. For implausibility reasons, we never choose initial infected people among persons in nursing homes or hospitals.

In each substep in each cycle contagion repeats its action, in different ways as we describe for the dirrefet sutuations and places.

We can set:

 * min and max duration of the infection;
 * the length of the incubation interval;
 * the critical distance, as the radius of a circle affecting people which are in it, with a given probability (how to set the value of the radius? See *Considerations on dimensions* below); 
 * the probability of the previous item, via *probabilityOfGettingInfection* corrected by *D%* > 0 if asymptomatic Infected people increase the probability of getting the infection because Susceptible subject go closer to them; *D%* < 0 if we suppose that asymptomatic Infected people are diffusing the infection in a lower way; *D%* = 0 means "no effects".

We have two types of contagion: (a) in a radius, for people moving around, also if only temporary present in a house/nursing home/hospital (in school we have only students and teacher). (b) in a given space (room or apartment) for people resident in their home /in a hospital / in a nursing home / or in being in school.

Resident people in hospitals and nursing homes suffer from two contagion possibilities: (a) and (b).

While people are at school, diffuse the infection to people in the same classroom, where only teachers and students are present, so this is a case (c).

### intrinsicalSusceptibilityFactor

We have also a subjective element set in *intrinsicalSusceptibilityFactor*. If *intrinsicalSusceptibilityFactor* is equal to 1, it has no effect. If > 1 it works, because: we will multiply the *probabilityOfGettingInfection*:

 * by 1 for fragile people; 
 * by *intrinsicalSusceptibilityFactor* for fragile people in nursing homes; 
 * by 1/*intrinsicalSusceptibilityFactor* for regular people;
 * by 1/*intrinsicalSusceptibilityFactor*<sup>2</sup> for robust people.

Never use *intrinsicalSusceptibilityFactor* < 1.

Both infected symptomatic (red) and infected asymptomatic (violet) go around infecting other people. 

Infected finally recover or die; their infection duration is supposed to be the same both for asymptomatic people and symptomatic ones.

Gray people with shape "circle" are regular persons; if they are young ones, 0-24 years, the internal variable *young* of the agents is set to True. Following data ([here](https://www.tuttitalia.it/piemonte/statistiche/popolazione-eta-sesso-stato-civile-2019/)) we set the number of young people to 21%.

Fragile people (gray with shape "x") have a high probability of becoming symptomatic, while regular people have a low probability of becoming symptomatic. In an area like that of the  Piedmont region, we have 25% pop. > 65 years (data [here](https://www.tuttitalia.it/piemonte/statistiche/popolazione-eta-sesso-stato-civile-2019/)). To take in consideration also the fragility among active workers with more than 55 years, we can increase the value at 30-35%.

Gray people with shape "star" or %triangle", are healthcare operators jumping in and out magenta and orange zone, so with high risk to be infected.

### City and Town effect

If we set the *cities&towns* switch to **On** (before Setup), the standard deviation of the movements of the persons in the white zone is reduced if a person is in a T or to a C patch (in the second case, more reduced).

The first effect is the emergence of zones more filled with people; the second one is a very different epidemic development, maybe more realistic.

Cities and towns also help us to simulate the existence of working areas.

### Apartments

To simulate the limitations/Lockdown, each house is sub-dived in 6 apartments. We have a number of houses close to 500. With 6 apartments in each house we have about 3,000 apartments; a number adequate to a population of 4,350 subjects. The actual number of apartments in Piedmont (4.35 millions of persons) is close to 2.8 millions (we can download a table [here](https://www.regione.piemonte.it/web/sites/default/files/media/documenti/2020-01/t13_15i.xls)). We set *populationSize* = 4350.

People at home stay in *their apartments*; they go home once each tick; then, at the beginning of the next tick, they to their previous positions.

### Population, houses, nursing homes, and hospitals


*Population*

All agents have their home; the procedure *createPeople* assigns a home to each agent, inside a city or a town. The agents have also a regular place (RP) where they act and interact, moving around. These positions can be interpreted as free time elective places. Students and teachers, when we activate the school, have both RPs and the schools; healthcare operators have both RPs and hospital or nursing homes; finally, workers have both RPs and working places.

RP are mostly placed close to the home, inside the same city or town, or in a given small, medium or large radius. To rule people distribution the program uses a vector of radiuses, from where to extract a radius within which to place an agent.

Currently, *radiusVector* is [2 3 3 3 3 5 5 5 7 9], that means that in one case on ten the RP is in radius 2 from agent home, in four cases on ten in radius 3, etc.

About the movement of the agents, see below "Movements and position in each patch (square) of the world".

*Cities, towns and Houses*

Houses are always present in the world. We can activate or not nursing homes and hospitals. If *citie&towns* is off, idle hospitals and nursing homes are anyway present into the world, in lighter colors; if it is on, idle hospitals and nursing homes do not appear into the world.

With *cities&towns* **on**, the program sets 1.25% of the patches to be towns (T) and the 0.25% to be cities (C), so we have about the 1.5% of 51•51 (2601) => about 32 patches, on a probabilistic basis, as downtown (1/6 as cities, so 5-6); around each center we add other patches in a small radius; globally, about 400 patches. We set 50% of those patches to contain houses. The construction is based on probability and beside that the actual number of patches can be lower than the expected one, due to overlapping in construction among cities and towns; we can expect about 200 patches as houses. This value is close to the 8% of the patches. In Piedmont the surface of the houses is the 5.8% of the total, see in the Chapter 6.1 at p. 63 [here](https://www.regione.piemonte.it/web/sites/default/files/media/documenti/2019-02/monitoraggio_consumosuolo_2015.pdf).

With *cities&towns* **off**, the program directly constructs an area of houses in 4 rows of 51 patches, being 200/51 => 3.922 rounded to 4.

*Hospitals*

To have a reference point, the Piedmont with 4.35 milions of inhabitants has 47 [hospitals](https://www.regione.piemonte.it/web/temi/sanita/organizzazione-strutture-sanitarie/strutture-ospedaliere).

At our population scale, we should divide 47 by 1,000. To have a visible phenomenon, we create 5 *hospitals* in the left side of the first row of the world, with pink color (case *cities&towns* **off**) or 5 patches with pink color in the C positions, not overlapping the houses (case *cities&towns* **on**).


*Nursing homes*

To have a reference, the Piedmont - with 4.35 milions of inabitants - has 784 [nursing homes](http://www.regione.piemonte.it/cgi-bin/polsoc/ricerca/presidi/index.cgi), hosting the [3% of the population](https://tinyurl.com/y7er688n) with more than 65 years. The number of guests is about 33,000.

At our population scale, we should divide 784 by 1,000. To avoid an out of scale concentration of persons in a unique nursing home, we create 5 nursing home in the right side of the first row of the world, with orange color (case *cities&towns* **off**) or 5 patch with orange color positions in a T or C area (case *cities&towns* **on**).

If we have cities and towns, nursing homes are localized outside.

*nursing homes* will host the 3% of the fragile subjects of the model. In an area like that of the  Piedmont region, we have 25% pop. > 65 years (data [here](https://www.tuttitalia.it/piemonte/statistiche/popolazione-eta-sesso-stato-civile-2019/)). The number of guests would be about 3% * 25% * 4350 / 5 => about 6 per nursing home. We set *probOfBeingInNursingHome%* to 5 instead that 3 to have less subtle quantities of the phenomenon, with 5% * 25% * 4350 / 5 => about 10 persons per nursing home.

In each tick we can have or not (*peopleVisitingNHs* On or Off) people visiting the nursing homes. For visiting people the probability *visitingNHsP%* is set to 1.4%, giving about 65 persons, uniquely moving in the first of the daily movements. NHs have only been accessible in the first part of the epidemic and 65 visitors assure the presence of infected agents there, anyway not always but around in the 50% of the experiments. The presence/absence of visiting people in NHs is relevant to compare the outcomes of the two cases of experiments.

We also have nursing home operators and other persons working there with *numberOperatorsPerNH* set to 3. They always move, also in limitations/Lockdown periods.

### Schools

If *activateSchools* is On we see yellow patches representing schools. In Table 3 [here](https://www.miur.gov.it/documents/20182/0/Principali+dati+della+scuola+-+avvio+anno+scolastico+2019-2020.pdf/5c4e6cc5-5df1-7bb1-2131-884daf008088?version=1.0&t=1570015597058) we read that In Piedmont there are more than 3 thousand schools. Using the 1:1000 scale, we would have 3 schools in the simulation, but we are interested in the possibility of experimenting about school lock and unlock, so we decide to create 2 schools for each 5 towns and 1 for each city, with about 20 schools, each with one classroom or two (*howManyClass1oomsInASchool*). NB, we can set the number of classrooms uniquely at tick 0.

Schools are placed within urban spaces, in yellow patches where we build no houses. Only students and teachers can go there. Each school has 4 teachers. At the same address above, in Table 10 and 11 we read that the number of teachers is close to 60 thousand. We have about 50, so close in scale.

In Table 5 we have 526 thousand students, plus 100 thousand at the university level (other sources), so we have about the 14% of the Piedmont population in school, against the 21% of young people on population, so 2/3. Students will chose the school closest to their home and randomly the classroom.

600 students in 20 schools with 1 room, give the presence of about 30 per room; with 2 rooms, about 15.

Students and teachers do visit the nursing homes if *peopleVisitingNHs* allows it, being *on*.

*activateSchools* only works if *cities&Towns* is On.

If the activation at tick 0 is missing, it is impossible to activate schools in a successive tick.

*activateSchools* is stronger than limitations/Lockdown, students and theachers go to school in any case if °activateSchools* is on.

Students and teachers are placed randomly in the patch-school (in virtual classrooms, if any); the radius of infection of outside people, if very large, could interfere with them.

With %Students (%St in scripting), we limit the number of students going to school to the quota expressed here.

### Workplaces and workers (employees)

To create workplaces (factories, offices, shops), we use the data for Piedmont reported [here](http://www.unioncamere.gov.it/Atlante_2015/regioni/piemonte/index.html).

Synthetically: 

- productive units of any type and size of around 0.5 million;
- units with more than 50 employees are 2,700 with around 400,000 employees (140 per unit), 9% of the population (400,000/4,350,000);
- working people around 1.8 million, of which 1.4 million in small business units; Piedmont population: 4.35 million

In our artificial world, we build: 

- 5 "large" units, with 100 (top value 150) employees, give a total of 500 employees, slightly more than the 9% above;
- 100 small units with about ten employees (top value 15), giving a total of 1,000 employee

The actual number of the employees in each unit is randomly arising form initial setup: we create the worker status as an attribute of regular not young people and not teachers (about 1810) and of fragile people (about 1500 if we set 35 as *fragilePeople%*), choosing 86% of the first group (obtaining about 1560) and 16% of the second group (around 240), with 1800 as the total. Young people are supposed to be students or unemployed.

We create large units in a radius 10 from urban locations and small units in a radius 5. Workers choose a large unit with probability 1/3; else, they choose a small one; in both cases, the closer to their home.

Factories can be close or open following *%openFactoriesAlsoWhenLimitationsOn*; "open" allows a worker to move to the factory also in case of limitations/Lockdown.

If the activation at tick 0 is missing, it is impossible to activate factories in a successive tick. It is impossible deactivate activated factories, use *%openFactoriesAlsoWhenLimitationsOn*, in case also set to 0.

If *fragileWorkersAtHome* is *On* the total of the workers is unchanged, but the workers are all regular; the former fragile ones are replaced by regular workers. This is a counterfactual operation and we can activate it uniquely at tick 0. We suggest, in this case, of setting *Off* *peopleVisitingNHs*, to avoid that too much people of type fragile (not at work) are visiting the nursing homes, with an explosion of the epidemic.

The choice of *stopFragileWorkers* set to *On* works in another way. Fragile workers can move out of their homes following the *%PeopleAnyTypeNotSymptomaticLeavingHome* or the *%PeopleNotFragileNotSymptomaticLeavingHome* probabilities, but can go to work in no case. Remember that the regular case is that the workers (fragile or regular)  can go to their factory (if open) also when *limitations/Lockdown* choice is *On*.

## Absolute quarantine effect

We can obtain a separation effect related to fragile people in nursing homes if we stop both the activity of the operators (counterfactual) and people visiting the nursing homes.

We can create a quarantine for infected people if all are sent in hospitals, but without operators (counterfactual).

## Movements and position in each patch (square) of the world

Movements are random, with ∆x and ∆y components, extracting two values from a normal distribution to determine the arrival point. The standard deviation of the normal distribution is corrected with a *sdReductionF*, low in "T", high in "C"; in orange (nursig homes) and pink (hospitals) patches the reduction factor is that used for "T". 

When we send a person to a patch, the position will be at random around the center of the patch, but in case of the home (with a specific apartment) of the agent, the internal memory keeps track of the center of the patch, to identify other people there, having same home and apartment. In *nursing homes* there is a unique apartment.

Contagion is operating within *radiusOfInfection* distance, but at home or in a nursing home or hospital or school we consider all people who are resident there (same apartment/classromm/space).

### Movements of hospital and nursing home operators

Hospital and nursing home operators always move once per day, also in limitations/Lockdown periods. In the 70% of the cases they go to work, in the 30% they move around to any place.

### Parameters

*run#* is an instrumental variable used in BehaviorSpace to order the repetition of the runs; if the seed indicator is negative, *run#* is related to the actual seed used in that run (look at "Special values: < 0 and ≠ -9999" paragraph).

The crucial parameter is *populationSize*, that we set considering the size of the 
world where interactions and infections occur. 

At present we take as a reference the Piedmont region in Italy, with 4,350 thousands of inhabitants (data [here](https://www.tuttitalia.it/piemonte/statistiche/popolazione-eta-sesso-stato-civile-2019/)).

The frequency of actions in each tick, *#movementCyclesPerTick*, is strictly related the size of the population in determining the spread of the virus. 

## Considerations on dimensions

Piedmont surface 2,539,699 ha ([source: Chapter 6.1 at p. 63](https://www.regione.piemonte.it/web/sites/default/files/media/documenti/2019-02/monitoraggio_consumosuolo_2015.pdf)). About 25,400 km2 or 25,400*10^6 m2.

Population 4350•10^3 units ([source](https://www.tuttitalia.it/piemonte/statistiche/popolazione-eta-sesso-stato-civile-2019/)).

Population / m2 => 4350 / 25,400•10^3 => 0.000171 (A),

Model surface 51•51 => 2601 patches (side of the square = 1); subdividing each square in 9 equal smaller squares with side 1/3 => 0.333É and surface (1/3)^2 => 0.111É we have 2601•9 => 23408 small squares (sSqs) with the same total surface.

Model population 4350, so population / sSqs => 4350 / 23408 => 0.1858, which is close to be on scale with 1000 • A => 0.171 per space unit (B) of 1000 m2.

Units B are square with side 1000^0.5 => 31,62.

In the real world, the three subdivisions of a model, patch each of 0.333 units, would be of 31.62 meters each, so 1 as the side of a patch is equivalent to 3•31,62; let us approximate to 100 meters.

A radius of 0.2 (the current one) is 20 meters, so we have to better calibrate with movement and probability, this is the next step. 


## HOW TO USE IT

Familiarize with the effects of the parameters, to prepare to experiment with different choices. Examples: 

 • limitations/Lockdown activation / deactivation;

 • % of infected people sent to hospitals;

 • Protecton of fragile people in nursing homes, stopping visitors;

 • Different rates of people circulation in limitations/Lockdown and when releasing it.

## USING BEHAVIOR SPACE

(For advanced users.)

The program contains a set of instructions to use the tool *Behavior Space*, with names *experiments 2 seeds*, with 1000 repetition.

Pay attention to the definition of the seeds (see above)..

In *experiments 2 seeds* we have:

- in the first box "Vary vairiables as follows" the row
["run#" [1 1 1000]]
determining 1000 runs with run#, instrumental variable created just for this use, from 1 to 1000;

- in Repetitions we have 1, as the repetions are following run#;

- in the second box "Measure runs using these reporters" (condering the report of all the possible checkpoints):

	seed0
seed
cp1cRed
cp1cRedNH
cp1cViolet
cp1cVioletNH
cp1cDead
cp2cRed
cp2cRedNH
cp2cViolet
cp2cVioletNH
cp2cDead
cp3cRed
cp3cRedNH
cp3cViolet
cp3cVioletNH
cp3cDead
cp4cRed
cp4cRedNH
cp4cViolet
cp4cVioletNH
cp4cDead
cp5cRed
cp5cRedNH
cp5cViolet
cp5cVioletNH
cp5cDead
cp6cRed
cp6cRedNH
cp6cViolet
cp6cVioletNH
cp6cDead
cumulativeRedCount
cumulativeVioletCount
cumulativeDeadCount
ticks
TinfSch
SinfSch
TinfNotSch
SinfNotSch

- In "Stop condition":

	theEnd  (a logic *true/false* variable).


We suggest to download the program to observe the repetition of the experiments, with different starting points (seeds) of the random numbers. This operation is quite smooth, via Behavior Search, which is in the menu Tools of NetLogo. It is possible to analyze the different results of the model at the various checkpoints. 

***Nota bene***. In each series of experiments, we will have variability in the results, although internally consistent. The model makes many agents act and interact; from the [complexity](https://terna.to.it/simul/complexity.html) point of view of each execution produces a "story," with a specific and unique sequence of emerging effects. Real events would behave in the same way if it were possible to repeat them.


## Rt CALCULATION (uniquely for reporting purposes)

Instead of *R* with subscript t, we write here Rt. NB, Rt instead of R0, which has a meaning only in a static view, at the beginning. This is a measure on move.

Let us suppose, for simplicity, u = 0 and d = 0 and that the I are of a unique symptomatic type

N = S + I + R

We have s(t) = S(t)/N, with each infected individual having a fixed number of contact (with S(t) or I(t) people) **b** in each day, with fixed probabilty of passing the infection on the fraction s(t) of the contacted persons. Into the simulation, the prob. can be modified and it is incorporated in the calculations made for evaluating Rt.

b s(t) is the number of new infected per each infected individual (in the code of the model, newRedCount/cumulativeRedCount_t-1 ); we can also note that b I(t) s(t) is equivalent to newRedCount and that the action of the *probabilityOfGettingInfectioné is
incorporated in newRedCount

**k** is the fraction of I(t) recovering in each day (in the code of the model, newGreenCount/cumulativeRedCount_t-1 ). As an example, I(t) is represented in the calculations by cumulativeRedCount_t-1.

Rt = b/k

We have [NB, _t-1 means at time -1]

1) Rt for symptomatic cases, reported in red

            newRedCount
        ----------------------
        cumulativeRedCount_t-1      newRedCount
    Rt= ---------------------- =  -------------- ; the * newGreenCount* means "coming
           newGreenCount*         newGreenCount*   from red subjects uniquely"
        ----------------------
        cumulativeRedCount_t-1

2) Rt for symptomatic + asymptomatic cases, reported in brown

                   newRedCount + newVioletCount
        --------------------------------------------------     
        cumulativeRedCount_t-1 + cumulativeVioletCount_t-1    newRedCount+newVioletCount
    Rt= -------------------------------------------------- = ---------------------------
                        newGreenCount                              newGreenCount
        --------------------------------------------------
        cumulativeRedCount_t-1 + cumulativeVioletCount_t-1

The first indicator is employed in real world; the second one, more meaningful, cannot be calculated in real world. 

## SAVING THE CONTAGION SEQUENCE TO A CSV FILE

We activate the saving of the contagion sequence, adding '0 csv 1' in the script. The program will ask the name of the file; avoid to include the file extensionthe into the reply.

Each record of the file contains:

- the ID of the agent transmitting the contagion (for the initial cases, externally generated, the ID has value -1);
    
- its contagion progressive number, starting from 1 (for the initial cases, externally generated, this value is 0);

- the conventional color of the place where it turned infected, following the NetLogo color swatches (for the externally generated initial cases, this value is 0);

- the ID of the agent receiving the contagion;

- its fragility rate (1 - robust; 2 - regular; 3 - fragile; 4 - extra fragile;

- its progressive contagion number;

- the conventional color of the place where it is turning infected, following the NetLogo color swatches;

- the day (tick) of the contagion (for the initial cases, externally generated, this value is 0);

- the starting infection day, i.e., the previous value plus the *incubationPeriod* (for the externally generated initial cases, the starting infection value is 0);

- the day of the conclusion of the infection, i.e., the previous value plus a value between the  *minInfectionDuration* and the *maxInfectionDuration* settings; this period stops if the agent deceases, but we do not consider that possibility here;

- the symptomatic (1) or asymptomatic (2) status.

## VACCINATION MANAGER

Is the script containa 

    0 vacc name

where name is the name of a file "name.csv", the program activates the vaccination procedures.

In name.csv the initial two values are:

numberOfDaysForImmunityActivation probabilityForTheVaccinatedPersonToBeContagious
with the second with base 1.

### Using Genetic Algorithms

The basis of our work is the agent-based simulation of an epidemic with propagation generated by highly mutable individual agent contacts. This is a model that inherently generates high variability in epidemic trends because even rare sequences of contagions can lead to very different overall outcomes. For this reason, in comparative applications the simulation is performed in repetition batches of one thousand or ten thousand times and those considered are mean values and distributions.

The genetic algorithms algorithm uses only one of those cases, carefully chosen by ourselves as a representative case, but the modifications induced by the parameter changes caused by the GA can create outbreak cases that are very rare to occur, and also counterintuitive. An example is the cases of contagion that appear to increase in the presence of vaccination, but actually as a result of different internal sequences of events due to changes in the agents (e.g., an infected person not leaving quarantine or vice versa).

This problem luckily finds compensation inside of the execution of the GA, because the simulation model is executed hundreds or thousands of times, being the algorithm evolving a population of models. Within those repetitions the extreme cases compensate themselves and also they change quickly.

To activate the GAs mechanism put into the script a row

0 vacc "some_name.csv" 
qith a content as
"y",40,1
373,5,0.10,0.10,0.10,0.10,0.10,0.10,0.10
433,10,0.10,0.10,0.10,0.10,0.10,0.10,0.10
493,10,0.10,0.10,0.10,0.10,0.10,0.10,0.10
553,10,0.10,0.10,0.10,0.10,0.10,0.10,0.10
613,20,0.10,0.10,0.10,0.10,0.10,0.10,0.10
738,"end"

373 means from day 373 using 5 (1000) doses vaccinate the 10% of the people in the first column, etc. Explanantions at https://terna.to.it/simul/GAresultPresentation.pdf

GAs use the sliders on the far right of the interface, if "y" is the first item in the .csv file named in the vacc item into the script; in that case, also ourselves can use drectly those sliders to set vaccination quotas. The second item measures the lenght (in days) of the interval for the effectiveness of the vaccination; the third one is the probability of contagion transmission from a vaccinated person, from 0 to 1.

**Hard finish** button, key "h" (decide if useful]

## How to cite

Terna P., Pescarmona G., Acquadro A., Pescarmona P., Russo G., Terna S. (2020), An Agent-Based Model of the Diffusion of Covid-19 Using NetLogo, URL https://terna.to.it/simul/SIsaR.html

Corresponding author Pietro Terna [pietro.terna@unito.it](mailto:pietro.terna@o
unito.it)

## Acknowledgments.

Many thanks to Simone Landini, Nizar Mansour, Fabio Pammolli, Enrico Scalas, and Federico Tedeschi for precious discussions, insights, and critics. The usual disclaimer applies.


## COPYRIGHT AND LICENSE

Copyright Pietro Terna, Gianpiero Pescarmona, Alberto Acquadro, Paolo Pescarmona, Giuseppe Russo, Stefano Terna.


![CC BY-NC-SA 3.0](http://ccl.northwestern.edu/images/creativecommons/byncsa.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License.  To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

building institution
false
0
Rectangle -7500403 true true 0 60 300 270
Rectangle -16777216 true false 130 196 168 256
Rectangle -16777216 false false 0 255 300 270
Polygon -7500403 true true 0 60 150 15 300 60
Polygon -16777216 false false 0 60 150 15 300 60
Circle -1 true false 135 26 30
Circle -16777216 false false 135 25 30
Rectangle -16777216 false false 0 60 300 75
Rectangle -16777216 false false 218 75 255 90
Rectangle -16777216 false false 218 240 255 255
Rectangle -16777216 false false 224 90 249 240
Rectangle -16777216 false false 45 75 82 90
Rectangle -16777216 false false 45 240 82 255
Rectangle -16777216 false false 51 90 76 240
Rectangle -16777216 false false 90 240 127 255
Rectangle -16777216 false false 90 75 127 90
Rectangle -16777216 false false 96 90 121 240
Rectangle -16777216 false false 179 90 204 240
Rectangle -16777216 false false 173 75 210 90
Rectangle -16777216 false false 173 240 210 255
Rectangle -16777216 false false 269 90 294 240
Rectangle -16777216 false false 263 75 300 90
Rectangle -16777216 false false 263 240 300 255
Rectangle -16777216 false false 0 240 37 255
Rectangle -16777216 false false 6 90 31 240
Rectangle -16777216 false false 0 75 37 90
Line -16777216 false 112 260 184 260
Line -16777216 false 105 265 196 265

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

factory
false
0
Rectangle -7500403 true true 76 194 285 270
Rectangle -7500403 true true 36 95 59 231
Rectangle -16777216 true false 90 210 270 240
Line -7500403 true 90 195 90 255
Line -7500403 true 120 195 120 255
Line -7500403 true 150 195 150 240
Line -7500403 true 180 195 180 255
Line -7500403 true 210 210 210 240
Line -7500403 true 240 210 240 240
Line -7500403 true 90 225 270 225
Circle -1 true false 37 73 32
Circle -1 true false 55 38 54
Circle -1 true false 96 21 42
Circle -1 true false 105 40 32
Circle -1 true false 129 19 42
Rectangle -7500403 true true 14 228 78 270

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

house colonial
false
0
Rectangle -7500403 true true 270 75 285 255
Rectangle -7500403 true true 45 135 270 255
Rectangle -16777216 true false 124 195 187 256
Rectangle -16777216 true false 60 195 105 240
Rectangle -16777216 true false 60 150 105 180
Rectangle -16777216 true false 210 150 255 180
Line -16777216 false 270 135 270 255
Polygon -7500403 true true 30 135 285 135 240 90 75 90
Line -16777216 false 30 135 285 135
Line -16777216 false 255 105 285 135
Line -7500403 true 154 195 154 255
Rectangle -16777216 true false 210 195 255 240
Rectangle -16777216 true false 135 150 180 180

house ranch
false
0
Rectangle -7500403 true true 270 120 285 255
Rectangle -7500403 true true 15 180 270 255
Polygon -7500403 true true 0 180 300 180 240 135 60 135 0 180
Rectangle -16777216 true false 120 195 180 255
Line -7500403 true 150 195 150 255
Rectangle -16777216 true false 45 195 105 240
Rectangle -16777216 true false 195 195 255 240
Line -7500403 true 75 195 75 240
Line -7500403 true 225 195 225 240
Line -16777216 false 270 180 270 255
Line -16777216 false 0 180 300 180

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiments 2 seeds" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>theEnd</exitCondition>
    <metric>seed0</metric>
    <metric>seed</metric>
    <metric>seedD</metric>
    <metric>cp1cRed</metric>
    <metric>cp1cRedNH</metric>
    <metric>cp1cViolet</metric>
    <metric>cp1cVioletNH</metric>
    <metric>cp1cDead</metric>
    <metric>cp2cRed</metric>
    <metric>cp2cRedNH</metric>
    <metric>cp2cViolet</metric>
    <metric>cp2cVioletNH</metric>
    <metric>cp2cDead</metric>
    <metric>cp3cRed</metric>
    <metric>cp3cRedNH</metric>
    <metric>cp3cViolet</metric>
    <metric>cp3cVioletNH</metric>
    <metric>cp3cDead</metric>
    <metric>cp4cRed</metric>
    <metric>cp4cRedNH</metric>
    <metric>cp4cViolet</metric>
    <metric>cp4cVioletNH</metric>
    <metric>cp4cDead</metric>
    <metric>cp5cRed</metric>
    <metric>cp5cRedNH</metric>
    <metric>cp5cViolet</metric>
    <metric>cp5cVioletNH</metric>
    <metric>cp5cDead</metric>
    <metric>cp6cRed</metric>
    <metric>cp6cRedNH</metric>
    <metric>cp6cViolet</metric>
    <metric>cp6cVioletNH</metric>
    <metric>cp6cDead</metric>
    <metric>cumulativeRedCount</metric>
    <metric>cumulativeVioletCount</metric>
    <metric>cumulativeDeadCount</metric>
    <metric>ticks</metric>
    <metric>TinfSch</metric>
    <metric>SinfSch</metric>
    <metric>TinfNotSch</metric>
    <metric>SinfNotSch</metric>
    <steppedValueSet variable="run#" first="1" step="1" last="1000"/>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
