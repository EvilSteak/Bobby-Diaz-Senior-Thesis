adc_data_path = "C:\\ti\\radar_data\\adc_data.bin"

if (ar1.ChanNAdcConfig(1, 1, 0, 1, 1, 1, 1, 2, 1, 0) == 0) then
    WriteToLog("ChanNAdcConfig Success\n", "green")
else
    WriteToLog("ChanNAdcConfig failure\n", "red")
end

if (ar1.LPModConfig(0, 1) == 0) then
    WriteToLog("LPModConfig Success\n", "green")
else
    WriteToLog("LPModConfig failure\n", "red")
end

if (ar1.RfInit() == 0) then
    WriteToLog("RfInit Success\n", "green")
else
    WriteToLog("RfInit failure\n", "red")
end

RSTD.Sleep(1000)

if (ar1.DataPathConfig(1, 1, 0) == 0) then
    WriteToLog("DataPathConfig Success\n", "green")
else
    WriteToLog("DataPathConfig failure\n", "red")
end

if (ar1.LvdsClkConfig(1, 1) == 0) then
    WriteToLog("LvdsClkConfig Success\n", "green")
else
    WriteToLog("LvdsClkConfig failure\n", "red")
end

if (ar1.LVDSLaneConfig(0, 1, 1, 0, 0, 1, 0, 0) == 0) then
    WriteToLog("LVDSLaneConfig Success\n", "green")
else
    WriteToLog("LVDSLaneConfig failure\n", "red")
end

if (ar1.ProfileConfig(0, 77, 100, 6, 60, 0, 0, 0, 0, 0, 0, 29.982, 0, 256, 5000, 0, 0, 30) == 0) then
    WriteToLog("ProfileConfig Success\n", "green")
else
    WriteToLog("ProfileConfig failure\n", "red")
end

-- TDM chirp config
if (ar1.ChirpConfig(0, 0, 0, 0, 0, 0, 0, 1, 0, 0) == 0) then
    WriteToLog("ChirpConfig 0 (TX0) Success\n", "green")
else
    WriteToLog("ChirpConfig 0 failure\n", "red")
end

if (ar1.ChirpConfig(1, 1, 0, 0, 0, 0, 0, 0, 1, 0) == 0) then
    WriteToLog("ChirpConfig 1 (TX1) Success\n", "green")
else
    WriteToLog("ChirpConfig 1 failure\n", "red")
end

if (ar1.FrameConfig(0, 1, 8, 128, 50, 0, 0, 1) == 0) then
    WriteToLog("FrameConfig Success\n", "green")
else
    WriteToLog("FrameConfig failure\n", "red")
end

-- Disable test source AFTER FrameConfig to prevent automatic re-enable
RSTD.Sleep(200)
ar1.DisableTestSource(0)
RSTD.Sleep(200)
ar1.DisableTestSource(0)
WriteToLog("Test Source Disabled\n", "green")

-- DCA1000 setup
if (ar1.SelectCaptureDevice("DCA1000") == 0) then
    WriteToLog("SelectCaptureDevice Success\n", "green")
else
    WriteToLog("SelectCaptureDevice failure\n", "red")
end

if (ar1.CaptureCardConfig_EthInit("192.168.33.30", "192.168.33.180", "12:34:56:78:90:12", 4096, 4098) == 0) then
    WriteToLog("CaptureCardConfig_EthInit Success\n", "green")
else
    WriteToLog("CaptureCardConfig_EthInit failure\n", "red")
end

if (ar1.CaptureCardConfig_Mode(1, 2, 1, 2, 3, 30) == 0) then
    WriteToLog("CaptureCardConfig_Mode Success\n", "green")
else
    WriteToLog("CaptureCardConfig_Mode failure\n", "red")
end

if (ar1.CaptureCardConfig_PacketDelay(25) == 0) then
    WriteToLog("CaptureCardConfig_PacketDelay Success\n", "green")
else
    WriteToLog("CaptureCardConfig_PacketDelay failure\n", "red")
end

-- Capture
ar1.CaptureCardConfig_StartRecord(adc_data_path, 0)
RSTD.Sleep(1000)
ar1.StartFrame()
RSTD.Sleep(5000)

WriteToLog("Capture complete: " .. adc_data_path .. "\n", "green")