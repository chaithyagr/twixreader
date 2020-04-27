from __future__ import print_function   # for python 2.7 compatibility
import os
import re
import sys
import time
import numpy as np


def Twix(infile):
    """Function for reading siemens twix raw data files."""

    if isinstance(infile, str):
        # assume that complete path is given
        if infile[-4:].lower() != '.dat':
            infile += '.dat'   # adds filetype ending to file
    else:
        # filename not a string, so assume that it is the MeasID
        measID = infile
        infile = [f for f in os.listdir('.') if re.search(
            r'^meas_MID0*' + str(measID) + '.*\.dat$', f)]
        if len(infile) == 0:
            print('error: .dat file with measID', measID, 'not found')
            raise ValueError
        elif len(infile) > 1:
            print('multiple files with measID', measID,
                  'found, choosing first occurence')
        infile = infile[0]

    infile = os.path.realpath(infile)

    twix_list = list()

    fread = open(infile, 'rb')
    fread.seek(0, os.SEEK_END)
    fileSize = np.uint64(fread.tell())
    fread.seek(0, os.SEEK_SET)
    firstInt, secondInt = np.fromfile(fread, dtype=np.uint32, count=2)
    # lazy software version check (VB or VD?)
    if (firstInt < 10000) and (secondInt <= 64):
        version = 'vd'
        print('Software version: VD (!?)')
        measID, fileID = np.fromfile(fread, dtype=np.uint32, count=2)
        measOffset = list()
        measLength = list()
        # number of different scans in file stored in 2nd int (wip, only
        # one supported for now)
        NScans = secondInt
        for _ in range(NScans):
            offset, length = np.fromfile(fread, dtype=np.uint64, count=2)
            measOffset.append(offset)
            measLength.append(length)
            fread.seek(152 - 16, os.SEEK_CUR)
            twix_list.append(dict())
    else:
        # in VB versions, the first 4 bytes indicate the beginning of the
        # raw data part of the file
        version = 'vb'
        print('Software  : VB (!?)')
        # VB does not support multiple scans in one file:
        NScans = 1
        measOffset = [np.uint64(0)]
        measLength = [fileSize]
        twix_list.append(dict())

    for s in range(NScans):
        scanStart = measOffset[s]
        scanEnd = scanStart + measLength[s]
        cPos = measOffset[s]
        fread.seek(cPos, os.SEEK_SET)
        hdr_len = np.fromfile(fread, dtype=np.uint32, count=1)[0]
        # TODO: parse header information
        cPos = measOffset[s] + np.uint64(hdr_len)
        scanStart = cPos

        print('\nscan ', s)
        update_progress(cPos - scanStart, scanEnd - scanStart, True)
        while cPos + 128 < scanEnd:  # fail-safe not to miss MDH_ACQEND
            update_progress(cPos - scanStart, scanEnd - scanStart, False)
            fread.seek(cPos, os.SEEK_SET)

            mdh, flags = __evalMDH(fread, version)

            if mdh.ulDMALength == 0 or flags.MDH_ACQEND:
                cPos += mdh.ulDMALength
                # jump to next full 512 bytes
                if cPos % 512:
                    cPos += np.uint32(512 - cPos % 512)
                break
            elif flags.MDH_SYNCDATA:
                # skip SYNCDATA
                cPos += mdh.ulDMALength
                continue

            mdh_types = list()
            if flags.MDH_IMASCAN:
                mdh_types.append('ima')
            if flags.MDH_NOISEADJSCAN:
                mdh_types.append('noise')
            if flags.MDH_PHASCOR and not flags.MDH_PATREFSCAN:
                mdh_types.append('phasecorr')
            if flags.MDH_PHASESTABSCAN or flags.MDH_REFPHASESTABSCAN:
                mdh_types.append('phasestab')
            if (not flags.MDH_PHASCOR and flags.MDH_PATREFSCAN) or flags.MDH_PATREFANDIMASCAN:
                mdh_types.append('refscan')
            if flags.MDH_PATREFSCAN and flags.MDH_PHASCOR:
                mdh_types.append('refscanPC')
            if flags.MDH_RTFEEDBACK or flags.MDH_HPFEEDBACK:
                mdh_types.append('RTfeedback')

            for item in mdh_types:
                if item not in twix_list[s]:
                    twix_list[s][item] = generic_twix_object(infile, mdh)
                twix_list[s][item].checkMDH(mdh, flags, cPos)

            # jump to mdh of next scan
            cPos += mdh.ulDMALength

        for item in twix_list[s]:
            twix_list[s][item].clean(mdh)

    fread.close()

    return twix_list


def __evalMDH(fread, version='vd'):
    mdh = SimpleBunch()
    flags = SimpleBunch()

    if version == 'vb':
        mdh.version = 'vb'
        mdh.szScanHeader = np.uint64(0)  # [bytes]
        mdh.szChannelHeader = np.uint64(128)  # [bytes]
        blob = fread.read(mdh.szChannelHeader)
        bit_mask = blob[20]
    elif version == 'vd':  # vd/ve
        mdh.version = 'vd'
        # vd/ve introduces shorter headers for later channels
        mdh.szScanHeader = np.uint64(192)  # [bytes]
        mdh.szChannelHeader = np.uint64(32)  # [bytes]
        blob = fread.read(mdh.szScanHeader)
        bit_mask = blob[40]
    else:
        print('Error: Version', version, 'unknown')
        raise ValueError

    flagsAndDMALength = np.frombuffer(blob[:4], dtype=np.uint32)[0]
    mdh.ulDMALength = np.uint32(flagsAndDMALength % (2**25))
    flags.MDH_ACQEND = bool(bit_mask & 2**0)
    flags.MDH_SYNCDATA = bool(bit_mask & 2**5)

    if flags.MDH_SYNCDATA or flags.MDH_ACQEND or mdh.ulDMALength == 0:
        return mdh, flags

    mdh.ulPackBit = int(flagsAndDMALength & (2**26))
    mdh.ulPCI_rx = int(flagsAndDMALength // (2**26))
    mdh.measUID = np.frombuffer(blob[4:8], dtype=np.int32)[0]
    mdh.scanCounter = np.frombuffer(blob[8:12], dtype=np.uint32)[0]
    mdh.timeStamp = np.frombuffer(blob[12:16], dtype=np.uint32)[0]
    mdh.pmuTimeStamp = np.frombuffer(blob[16:20], dtype=np.uint32)[0]

    if mdh.version == 'vb':
        mdh.evalInfoMask = np.frombuffer(blob[20:28], dtype=np.uint32)
        mdh.samplesInScan = np.frombuffer(blob[28:30], dtype=np.uint16)[0]
        mdh.usedChannels = np.frombuffer(blob[30:32], dtype=np.uint16)[0]
        mdh.sLC = np.frombuffer(blob[32:60], dtype=np.uint16)
        mdh.cutOff = np.frombuffer(blob[60:64], dtype=np.uint16)
        mdh.kSpaceCentreColumn = np.frombuffer(blob[64:66], dtype=np.uint16)[0]
        mdh.coilSelect = np.frombuffer(blob[66:68], dtype=np.uint16)[0]
        mdh.readOutOffcentre = np.frombuffer(blob[68:72], dtype=np.float32)[0]
        mdh.timeSinceLastRF = np.frombuffer(blob[72:76], dtype=np.uint32)[0]
        mdh.kSpaceCentreLineNo = np.frombuffer(blob[76:78], dtype=np.uint16)[0]
        mdh.kSpaceCentrePartitionNo = np.frombuffer(blob[78:80], dtype=np.uint16)[0]
        mdh.iceProgramPara = np.frombuffer(blob[80:88], dtype=np.uint16)
        mdh.freePara = np.frombuffer(blob[88:96], dtype=np.uint16)
        mdh.sliceData = np.frombuffer(blob[96:124], dtype=np.float32)
        mdh.channelID = np.frombuffer(blob[124:126], dtype=np.uint16)
        mdh.pTABPosNeg = np.frombuffer(blob[126:128], dtype=np.uint16)
    else:
        mdh.systemType = np.frombuffer(blob[20:22], dtype=np.uint16)[0]
        mdh.pTABPosDelay = np.frombuffer(blob[22:24], dtype=np.uint16)[0]
        mdh.pTABPosXYZ = np.frombuffer(blob[24:36], dtype=np.int32)
        mdh.Reserved1 = np.frombuffer(blob[36:40], dtype=np.uint32)[0]
        mdh.evalInfoMask = np.frombuffer(blob[40:48], dtype=np.uint32)
        mdh.samplesInScan = np.frombuffer(blob[48:50], dtype=np.uint16)[0]
        mdh.usedChannels = np.frombuffer(blob[50:52], dtype=np.uint16)[0]
        mdh.sLC = np.frombuffer(blob[52:80], dtype=np.uint16)
        mdh.cutOff = np.frombuffer(blob[80:84], dtype=np.uint16)
        mdh.kSpaceCentreColumn = np.frombuffer(blob[84:86], dtype=np.uint16)[0]
        mdh.coilSelect = np.frombuffer(blob[86:88], dtype=np.uint16)[0]
        mdh.readOutOffcentre = np.frombuffer(blob[88:92], dtype=np.float32)[0]
        mdh.timeSinceLastRF = np.frombuffer(blob[92:96], dtype=np.uint32)[0]
        mdh.kSpaceCentreLineNo = np.frombuffer(blob[96:98], dtype=np.uint16)[0]
        mdh.kSpaceCentrePartitionNo = np.frombuffer(blob[98:100], dtype=np.uint16)[0]
        mdh.sliceData = np.frombuffer(blob[100:128], dtype=np.float32)
        mdh.iceProgramPara = np.frombuffer(blob[128:176], dtype=np.uint16)
        # renamed reservedPara to freePara for compatibility to VB:
        mdh.freePara = np.frombuffer(blob[176:184], dtype=np.uint16)
        mdh.applicationCounter = np.frombuffer(blob[184:186], dtype=np.uint16)[0]
        mdh.applicationMask = np.frombuffer(blob[186:188], dtype=np.uint16)[0]
        mdh.CRC = np.frombuffer(blob[188:192], dtype=np.uint32)[0]

    # now evaluate sLC """
    mdh.lin = mdh.sLC[0]
    mdh.ave = mdh.sLC[1]
    mdh.sli = mdh.sLC[2]
    mdh.par = mdh.sLC[3]
    mdh.eco = mdh.sLC[4]
    mdh.phs = mdh.sLC[5]
    mdh.rep = mdh.sLC[6]
    mdh.set = mdh.sLC[7]
    mdh.seg = mdh.sLC[8]
    mdh.ida = mdh.sLC[9]
    mdh.idb = mdh.sLC[10]
    mdh.idc = mdh.sLC[11]
    mdh.idd = mdh.sLC[12]
    mdh.ide = mdh.sLC[13]

    # now evaluate evalInfoMask """
    # const MdhBitField MDH_ACQEND            ((unsigned long)0);
    flags.MDH_ACQEND = bool(mdh.evalInfoMask[0] & 2**0)
    # const MdhBitField MDH_RTFEEDBACK        (1);
    flags.MDH_RTFEEDBACK = bool(mdh.evalInfoMask[0] & 2**1)
    # const MdhBitField MDH_HPFEEDBACK        (2);
    flags.MDH_HPFEEDBACK = bool(mdh.evalInfoMask[0] & 2**2)
    # const MdhBitField MDH_ONLINE            (3);
    flags.MDH_ONLINE = bool(mdh.evalInfoMask[0] & 2**3)
    # const MdhBitField MDH_OFFLINE           (4);
    flags.MDH_OFFLINE = bool(mdh.evalInfoMask[0] & 2**4)
    # const MdhBitField MDH_SYNCDATA          (5);       // readout
    # contains synchroneous data
    flags.MDH_SYNCDATA = bool(mdh.evalInfoMask[0] & 2**5)
    # const MdhBitField MDH_LASTSCANINCONCAT  (8);       // Flag for last
    # scan in concatination
    flags.MDH_LASTSCANINCONCAT = bool(mdh.evalInfoMask[0] & 2**8)
    # const MdhBitField MDH_RAWDATACORRECTION (10);      // Correct the
    # rawadata with the rawdata correction factor
    flags.MDH_RAWDATACORRECTION = bool(mdh.evalInfoMask[0] & 2**10)
    # const MdhBitField MDH_LASTSCANINMEAS    (11);      // Flag for last
    # scan in measurement
    flags.MDH_LASTSCANINMEAS = bool(mdh.evalInfoMask[0] & 2**11)
    # const MdhBitField MDH_SCANSCALEFACTOR   (12);      // Flag for scan
    # specific additional scale factor
    flags.MDH_SCANSCALEFACTOR = bool(mdh.evalInfoMask[0] & 2**12)
    # const MdhBitField MDH_2NDHADAMARPULSE   (13);      // 2nd RF
    # exitation of HADAMAR
    flags.MDH_2NDHADAMARPULSE = bool(mdh.evalInfoMask[0] & 2**13)
    # const MdhBitField MDH_REFPHASESTABSCAN  (14);      // reference phase
    # stabilization scan
    flags.MDH_REFPHASESTABSCAN = bool(mdh.evalInfoMask[0] & 2**14)
    # const MdhBitField MDH_PHASESTABSCAN     (15);      // phase
    # stabilization scan
    flags.MDH_PHASESTABSCAN = bool(mdh.evalInfoMask[0] & 2**15)
    # const MdhBitField MDH_D3FFT             (16);      // execute 3D FFT
    flags.MDH_D3FFT = bool(mdh.evalInfoMask[0] & 2**16)
    # const MdhBitField MDH_SIGNREV           (17);      // sign reversal
    flags.MDH_SIGNREV = bool(mdh.evalInfoMask[0] & 2**17)
    # const MdhBitField MDH_PHASEFFT          (18);      // execute phase
    # fft
    flags.MDH_PHASEFFT = bool(mdh.evalInfoMask[0] & 2**18)
    # const MdhBitField MDH_SWAPPED           (19);      // swapped
    # phase/readout direction
    flags.MDH_SWAPPED = bool(mdh.evalInfoMask[0] & 2**19)
    # const MdhBitField MDH_POSTSHAREDLINE    (20);      // shared line
    flags.MDH_POSTSHAREDLINE = bool(mdh.evalInfoMask[0] & 2**20)
    # const MdhBitField MDH_PHASCOR           (21);      // phase
    # correction data
    flags.MDH_PHASCOR = bool(mdh.evalInfoMask[0] & 2**21)
    # const MdhBitField MDH_PATREFSCAN        (22);      // additonal scan
    # for PAT reference line/partition
    flags.MDH_PATREFSCAN = bool(mdh.evalInfoMask[0] & 2**22)
    # const MdhBitField MDH_PATREFANDIMASCAN  (23);      // additonal scan
    # for PAT reference line/partition that is also used as image scan
    flags.MDH_PATREFANDIMASCAN = bool(mdh.evalInfoMask[0] & 2**23)
    # const MdhBitField MDH_REFLECT           (24);      // reflect line
    flags.MDH_REFLECT = bool(mdh.evalInfoMask[0] & 2**24)
    # const MdhBitField MDH_NOISEADJSCAN      (25);      // noise adjust
    # scan --> Not used in NUM4
    flags.MDH_NOISEADJSCAN = bool(mdh.evalInfoMask[0] & 2**25)
    # const MdhBitField MDH_SHARENOW          (26);      // all lines are
    # acquired from the actual and previous e.g. phases
    flags.MDH_SHARENOW = bool(mdh.evalInfoMask[0] & 2**26)
    # const MdhBitField MDH_LASTMEASUREDLINE  (27);      // indicates that
    # the current line is the last measured line of all succeeding e.g.
    # phases
    flags.MDH_LASTMEASUREDLINE = bool(mdh.evalInfoMask[0] & 2**27)
    # const MdhBitField MDH_FIRSTSCANINSLICE  (28);      // indicates first
    # scan in slice (needed for time stamps)
    flags.MDH_FIRSTSCANINSLICE = bool(mdh.evalInfoMask[0] & 2**28)
    # const MdhBitField MDH_LASTSCANINSLICE   (29);      // indicates  last
    # scan in slice (needed for time stamps)
    flags.MDH_LASTSCANINSLICE = bool(mdh.evalInfoMask[0] & 2**29)
    # const MdhBitField MDH_TREFFECTIVEBEGIN  (30);      // indicates the
    # begin time stamp for TReff (triggered measurement)
    flags.MDH_TREFFECTIVEBEGIN = bool(mdh.evalInfoMask[0] & 2**30)
    # const MdhBitField MDH_TREFFECTIVEEND    (31);      // indicates the
    # end time stamp for TReff (triggered measurement)
    flags.MDH_TREFFECTIVEEND = bool(mdh.evalInfoMask[0] & 2**31)
    # we add another flag for imascans (just to simplify stuff later on)
    flags.MDH_IMASCAN = True

    if flags.MDH_RTFEEDBACK or flags.MDH_HPFEEDBACK or flags.MDH_PHASCOR or flags.MDH_NOISEADJSCAN:
        flags.MDH_IMASCAN = False

    # otherwise the PATREFSCAN might be overwritten
    if flags.MDH_PHASESTABSCAN or flags.MDH_REFPHASESTABSCAN:
        flags.MDH_PATREFSCAN = False
        flags.MDH_PATREFANDIMASCAN = False
        flags.MDH_IMASCAN = False

    if flags.MDH_PATREFSCAN and not flags.MDH_PATREFANDIMASCAN:
        flags.MDH_IMASCAN = False

    # pehses: the pack bit indicates that multiple ADC are packed into one
    # DMA, often in EPI scans (controlled by fRTSetReadoutPackaging in IDEA)
    # since this code assumes one adc (x NCha) per DMA, we have to correct
    # the "DMA length"
    # if mdh.ulPackBit:
    # it seems that the packbit is not always set correctly
    mdh.ulDMALength = np.uint32(mdh.szScanHeader +
                                (2 * 4 * mdh.samplesInScan + mdh.szChannelHeader) * mdh.usedChannels)

    return mdh, flags


class generic_twix_object(object):
    """generic data object - can be ima, refscan,..."""

    def __init__(self, infile, mdh):
        self.nCol = []
        self.nCha = []
        self.lin = []
        self.par = []
        self.sli = []
        self.ave = []
        self.phs = []
        self.eco = []
        self.rep = []
        self.set = []
        self.seg = []
        self.ida = []
        self.idb = []
        self.idc = []
        self.idd = []
        self.ide = []
        self.centCol = []
        self.centLin = []
        self.centPar = []
        self.isReflected = []
        self.needsRawDataCorr = []
        self.memPos = []
        self.sliceData = []
        self.icePara = []
        self.freePara = []
        self.nAcq = 0
        self.infile = infile
        # save some info from mdh
        self.version = mdh.version
        self.szScanHeader = mdh.szScanHeader
        self.szChannelHeader = mdh.szChannelHeader
        self.dataDims = ['Col', 'Cha', 'Lin', 'Par', 'Sli', 'Ave',
                         'Phs', 'Eco', 'Rep', 'Set', 'Seg', 'Ida',
                         'Idb', 'Idc', 'Idd', 'Ide']
        # self.args = SimpleBunch()
        self.args = SimpleBunchNotify(self)
        self.args.averageFlags = {d: False for d in self.dataDims}
        self.args.removeOS = False
        self.args.skipLeadingZeros = False
        # options: colFirst, colLast, colLast is more numpy-style
        self.args.dataOrder = 'colLast'

    def checkMDH(self, mdh, flags, cPos):
        # save mdh information about current line
        self.nAcq += 1
        self.nCol.append(mdh.samplesInScan)
        self.nCha.append(mdh.usedChannels)
        self.lin.append(mdh.lin)
        self.par.append(mdh.par)
        self.sli.append(mdh.sli)
        self.ave.append(mdh.ave)
        self.phs.append(mdh.phs)
        self.eco.append(mdh.eco)
        self.rep.append(mdh.rep)
        self.set.append(mdh.set)
        self.seg.append(mdh.seg)
        self.ida.append(mdh.ida)
        self.idb.append(mdh.idb)
        self.idc.append(mdh.idc)
        self.idd.append(mdh.idd)
        self.ide.append(mdh.ide)
        self.centCol.append(mdh.kSpaceCentreColumn)
        self.centLin.append(mdh.kSpaceCentreLineNo)
        self.centPar.append(mdh.kSpaceCentrePartitionNo)
        self.isReflected.append(bool(flags.MDH_REFLECT))
        self.needsRawDataCorr.append(bool(flags.MDH_RAWDATACORRECTION))
        self.sliceData.append(mdh.sliceData)
        self.icePara.append(mdh.iceProgramPara)
        self.freePara.append(mdh.freePara)
        # save memory position
        self.memPos.append(cPos)

    def clean(self, mdh):
        # cast lists to arrays
        self.lin = np.asarray(self.lin, np.uint32)
        self.par = np.asarray(self.par, np.uint32)
        self.sli = np.asarray(self.sli, np.uint32)
        self.ave = np.asarray(self.ave, np.uint32)
        self.phs = np.asarray(self.phs, np.uint32)
        self.eco = np.asarray(self.eco, np.uint32)
        self.rep = np.asarray(self.rep, np.uint32)
        self.set = np.asarray(self.set, np.uint32)
        self.seg = np.asarray(self.seg, np.uint32)
        self.ida = np.asarray(self.ida, np.uint32)
        self.idb = np.asarray(self.idb, np.uint32)
        self.idc = np.asarray(self.idc, np.uint32)
        self.idd = np.asarray(self.idd, np.uint32)
        self.ide = np.asarray(self.ide, np.uint32)
        self.nLin = max(self.lin) + 1
        self.nPar = max(self.par) + 1
        self.nSli = max(self.sli) + 1
        self.nAve = max(self.ave) + 1
        self.nPhs = max(self.phs) + 1
        self.nEco = max(self.eco) + 1
        self.nRep = max(self.rep) + 1
        self.nSet = max(self.set) + 1
        self.nSeg = max(self.seg) + 1
        self.nIda = max(self.ida) + 1
        self.nIdb = max(self.idb) + 1
        self.nIdc = max(self.idc) + 1
        self.nIdd = max(self.idd) + 1
        self.nIde = max(self.ide) + 1
        self.memPos = np.asarray(self.memPos)

        # size of one "dataset" (i.e nCha*nCol + nCha*szChannelheader)
        self.readSz = [self.nCha[0],
                       (int(self.szChannelHeader // 8) + self.nCol[0])]
        self.__updateIndices()

    def notifyChange(self):
        self.__updateIndices()

    def __updateIndices(self):
        self.dataDims = ['Col', 'Cha', 'Lin', 'Par', 'Sli', 'Ave',
                         'Phs', 'Eco', 'Rep', 'Set', 'Seg', 'Ida',
                         'Idb', 'Idc', 'Idd', 'Ide']

        self.fullSize = [self.nCol[0], self.nCha[0], self.nLin, self.nPar,
                         self.nSli, self.nAve, self.nPhs, self.nEco,
                         self.nRep, self.nSet, self.nSeg, self.nIda,
                         self.nIdb, self.nIdc, self.nIdd, self.nIde]

        if self.args.dataOrder.lower() == 'collast':
            self.dataDims.reverse()
            self.fullSize = self.fullSize[::-1]

        self.shape = self.fullSize.copy()

        if self.args.removeOS:
            self.shape[self.dataDims.index('Col')] = self.nCol[0] // 2

        if self.args.skipLeadingZeros:
            self.skipLin = min(self.lin)
            self.skipPar = min(self.par)
        else:
            self.skipLin = 0
            self.skipPar = 0

        self.shape[self.dataDims.index('Lin')] = self.nLin - self.skipLin
        self.shape[self.dataDims.index('Par')] = self.nPar - self.skipPar

        # calculate indices to target & source(raw)
        linIx = self.lin - self.skipLin
        parIx = self.par - self.skipPar
        nLinIx = self.nLin - self.skipLin
        nParIx = self.nPar - self.skipPar

        if self.args.dataOrder.lower() == 'colfirst':
            self.ixToTarget = self.ide + self.nIde * (self.idd + self.nIdd * (self.idc + self.nIdc *
                                                                              (self.idb + self.nIdb * (self.ida + self.nIda * (self.seg + self.nSeg *
                                                                                                                               (self.set + self.nSet * (self.rep + self.nRep * (self.eco + self.nEco *
                                                                                                                                                                                (self.phs + self.nPhs * (self.ave + self.nAve * (self.sli + self.nSli *
                                                                                                                                                                                                                                 (parIx + nParIx * (linIx)))))))))))))
            # now calculate inverse index
            # inverse indices of lines that are not measured are -1
            self.ixToRaw = -np.ones(np.prod(self.shape[2:]), np.int32)
            for k in range(len(self.ixToTarget)):
                self.ixToRaw[self.ixToTarget[k]] = k
            self.ixToRaw = self.ixToRaw.reshape(self.shape[2:])
        elif self.args.dataOrder.lower() == 'collast':
            self.ixToTarget = linIx + nLinIx * (parIx + nParIx * (self.sli + self.nSli *
                                                                  (self.ave + self.nAve * (self.phs + self.nPhs * (self.eco + self.nEco *
                                                                                                                   self.rep + self.nRep * (self.set + self.nSet * (self.seg + self.nSeg *
                                                                                                                                                                   self.ida + self.nIda * (self.idb + self.nIdb * (self.idc + self.nIdc *
                                                                                                                                                                                                                   self.idd + self.nIdd * (self.ide))))))))))
            # now calculate inverse index
            # inverse indices of lines that are not measured are -1
            self.ixToRaw = -np.ones(np.prod(self.shape[:-2]), np.int32)
            for k in range(len(self.ixToTarget)):
                self.ixToRaw[self.ixToTarget[k]] = k
            self.ixToRaw = self.ixToRaw.reshape(self.shape[:-2])
        else:
            raise ValueError

    def __getitem__(self, index):
        """Overloading of the [] operator."""
        # translate index for array of size self.shape
        # to index for array of size self.shape[:-2]
        ndims = len(self.shape)
        nindex = np.size(index)
        if nindex == 1:
            index = [index]
        newIndex = list()
        for count in range(nindex):
            ix = index[count]
            remaining = nindex - count - 1
            if ix == Ellipsis:
                newIndex.extend([slice(None, None, None)] *
                                (ndims - len(newIndex) - remaining))
            else:
                newIndex.append(ix)
        # now fill remaining dims with ':'
        newIndex.extend([slice(None, None, None)] * (ndims - len(newIndex)))
        index = newIndex

        selCol = np.arange(self.shape[self.dataDims.index('Col')])[
            index[self.dataDims.index('Col')]]
        nSelCol = np.size(selCol)
        selCha = np.arange(self.shape[self.dataDims.index('Cha')])[
            index[self.dataDims.index('Cha')]]
        nSelCha = np.size(selCha)
        # select ixToRaw
        if self.args.dataOrder.lower() == 'colfirst':
            cIxToRaw = self.ixToRaw[index[2:]]
            outShape = list([nSelCol, nSelCha])
            outShape.extend(cIxToRaw.shape)
        elif self.args.dataOrder.lower() == 'collast':
            cIxToRaw = self.ixToRaw[index[:-2]]
            outShape = list(cIxToRaw.shape)
            outShape.extend([nSelCha, nSelCol])
        else:
            raise ValueError

        # vectorize ixToRaw
        cIxToRaw = cIxToRaw.flatten()

        return self.__read_data(cIxToRaw, outShape, index)

    def raw(self, sel=slice(None)):
        ncol = self.shape[self.dataDims.index('Col')]
        ncha = self.shape[self.dataDims.index('Cha')]
        nsel = len(np.zeros(len(self.memPos))[sel])
        if self.args.dataOrder.lower() == 'colfirst':
            outShape = list([ncol, ncha, nsel])
        elif self.args.dataOrder.lower() == 'collast':
            outShape = list([nsel, ncha, ncol])
        else:
            raise ValueError
        return self.__read_data(sel, outShape)

    def __read_data(self, cIxToRaw, outShape, index=[slice(None), slice(None)]):
        out = np.zeros(outShape, np.complex64)
        if self.args.dataOrder.lower() == 'colfirst':
            out = np.reshape(out, (outShape[0], outShape[1], -1))
        elif self.args.dataOrder.lower() == 'collast':
            out = np.reshape(out, [-1, outShape[-2], outShape[-1]])
        else:
            raise ValueError

        if self.args.removeOS:
            cutOS = np.concatenate([np.arange(self.nCol[0] // 4),
                                    np.arange(self.nCol[0] * 3 // 4, self.nCol[0])])

        if type(cIxToRaw) is slice:
            raw_mode = True
            # select memory
            mem = self.memPos[cIxToRaw].copy()
        else:
            raw_mode = False
            # create inverted indices
            cIxToTarg = np.arange(len(cIxToRaw))
            # delete all '-1's in cIxToRaw (not acquired)
            acquired = np.where(cIxToRaw >= 0)
            cIxToRaw = cIxToRaw[acquired]
            cIxToTarg = cIxToTarg[acquired]

            # select memory
            mem = self.memPos[cIxToRaw].copy()

            # sort mem for faster access
            ix = np.argsort(mem)
            mem = mem[ix]
            cIxToRaw = cIxToRaw[ix]
            cIxToTarg = cIxToTarg[ix]

        fread = open(self.infile, 'rb')
        for k, pos in enumerate(mem):
            # jump to current mem. pos.
            fread.seek(pos + self.szScanHeader)

            # read data
            raw = np.fromfile(fread, dtype=np.complex64,
                              count=np.prod(self.readSz))

            raw = np.reshape(raw, self.readSz)
            # now cut channelheader
            raw = raw[:, int(self.szChannelHeader // 8):]

            # select channels
            raw = np.atleast_2d(raw[index[-2], :])

            if self.args.removeOS:
                raw = np.fft.fft(np.fft.ifft(raw)[:, cutOS])

            if raw_mode:
                trg_ix = k
                raw_ix = k
            else:
                trg_ix = cIxToTarg[k]
                raw_ix = cIxToRaw[k]

            if self.isReflected[raw_ix]:
                raw = raw[:, -1:: -1]

                # select col and sort data
            if self.args.dataOrder.lower() == 'colfirst':
                out[:, :, trg_ix] += raw[:, index[-1]].transpose()
            elif self.args.dataOrder.lower() == 'collast':
                out[trg_ix, :, :] += raw[:, index[-1]]
            else:
                raise ValueError

        fread.close()

        return np.reshape(out, outShape)


def update_progress(cPos, endPos, initialize=False):
    if initialize:
        update_progress.tstart = time.time()
        update_progress.preprogress = -1

    progress = round(float(cPos) / endPos * 100)
    if progress - update_progress.preprogress >= 1:
        elapsed_time = (time.time() - update_progress.tstart)
        remain_time = elapsed_time * (100. / max(progress, 1) - 1.)
        # sys.stdout.mode
        sys.stdout.write('\r%3d %% parsed in %d s. Approximately %d s remaining.' % (
            progress, round(elapsed_time), round(remain_time)))
        sys.stdout.flush()
        update_progress.preprogress = progress


class SimpleBunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class SimpleBunchNotify(object):
    """
    Simple class that holds an argument dictionary.

    Triggers a call to .notifyChange() of the observer class(es) whenever
    a key is changed (not when it is initialized)
    """

    def __init__(self, *observers, **kwds):
        self.__dict__.update(kwds)
        self.observes = observers

    def __setattr__(self, key, value):
        if key in self.__dict__:
            object.__setattr__(self, key, value)
            for obs in self.observes:
                obs.notifyChange()
        else:
            object.__setattr__(self, key, value)


def twixprot(infile):
    """
    Parses xprotocol of siemens twix-file (.dat) and returns the protocol as a SimpleBunch().

    Currently only (the first occurence of) the ascconv-protocol
    data is read. Seems to be enough information for now.
    """
    if isinstance(infile, str):
        # assume that complete path is given
        if infile[-4:].lower() != '.dat':
            infile += '.dat'   # adds filetype ending to file
    else:
        # filename not a string, so assume that it is the MeasID
        measID = infile
        infile = [f for f in os.listdir('.') if re.search(
            r'^meas_MID0*' + str(measID) + '.*\.dat$', f)]
        if len(infile) == 0:
            print('error: .dat file with measID', measID, 'not found')
            raise ValueError
        elif len(infile) > 1:
            print('multiple files with measID', measID,
                  'found, choosing first occurence')
        infile = infile[0]

    f = open(infile, 'rb')

    # find begin of assconv data
    while str(f.readline()).find('### ASCCONV BEGIN') == -1:
        pass

    ascconv = SimpleBunch()
    # ascconv = dict()
    for line in f.readlines():
        line = line.decode('utf-8').rstrip()
        if line.find("### ASCCONV END") > -1:
            break
        elif line:
            cBunch = ascconv
            parts = line.split()
            splitparts = parts[0].split('.')
            for i, var in enumerate(splitparts):
                if var.startswith('a'):
                    tmp = var.strip(']').split('[')
                    key = tmp[0]
                    if len(tmp) > 1:
                        index = int(tmp[1])
                    else:
                        index = 0
                    if key not in cBunch.keys():
                        cBunch.update({key: list()})
                    if i < len(splitparts) - 1:
                        for k in range(len(cBunch[key]) - 1, index):
                            cBunch.get(key).append(SimpleBunch())
                        cBunch = cBunch[key][index]
                    else:
                        val = parts[-1]
                        if key.startswith(('al', 'an', 'ac')):
                            val = int(float(val))
                        elif key.startswith(('ad', 'afl')):
                            val = float(val)
                        elif key.startswith(('auc', 'aui', 'aul', 'aun')):
                            val = int(val, 16)
                        else:
                            val = val.strip('"')
                        for k in range(len(cBunch[key]) - 1, index - 1):
                            cBunch.get(key).append([])
                        cBunch.get(key).insert(index, val)
                else:
                    if i < len(splitparts) - 1:
                        if var not in cBunch.keys():
                            cBunch.update({var: SimpleBunch()})
                        cBunch = cBunch[var]
                    else:
                        key = splitparts[-1]
                        val = parts[-1]
                        if key.startswith(('l', 'n', 'c')):
                            val = int(float(val))
                        elif key.startswith(('d', 'fl')):
                            val = float(val)
                        elif key.startswith(('uc', 'ui', 'ul', 'un', 'e', 'size')):
                            val = int(val, 16)
                        else:
                            val = val.strip('"')
                        cBunch.update({key: val})
    return ascconv
