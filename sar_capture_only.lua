PC_IP      = "192.168.33.30"
DCA_IP     = "192.168.33.180"
DCA_MAC    = "12:34:56:78:90:12"
CMD_PORT   = 4096
DATA_PORT  = 4098

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

-- TDM chirp config: chirp 0 = TX0 only, chirp 1 = TX1 only
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

-- 65535 frames at 50ms = ~54 minutes continuous transmission
if (ar1.FrameConfig(0, 1, 65535, 128, 50, 0, 0, 1) == 0) then
    WriteToLog("FrameConfig Success (65535 frames)\n", "green")
else
    WriteToLog("FrameConfig failure\n", "red")
end

-- Disable test source AFTER FrameConfig
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

if (ar1.CaptureCardConfig_EthInit(PC_IP, DCA_IP, DCA_MAC, CMD_PORT, DATA_PORT) == 0) then
    WriteToLog("CaptureCardConfig_EthInit Success\n", "green")
else
    WriteToLog("CaptureCardConfig_EthInit failure\n", "red")
end

if (ar1.CaptureCardConfig_Mode(1, 2, 1, 2, 3, 0) == 0) then
    WriteToLog("CaptureCardConfig_Mode Success (unlimited UDP stream)\n", "green")
else
    WriteToLog("CaptureCardConfig_Mode failure\n", "red")
end

if (ar1.CaptureCardConfig_PacketDelay(25) == 0) then
    WriteToLog("CaptureCardConfig_PacketDelay Success\n", "green")
else
    WriteToLog("CaptureCardConfig_PacketDelay failure\n", "red")
end

-- Start recording to trigger UDP stream
-- The file will grow but UDP packets flow simultaneously
-- Set a dummy path — we only care about the UDP stream in Python
if (ar1.CaptureCardConfig_StartRecord("C:\\ti\\radar_data\\stream.bin", 1) == 0) then
    WriteToLog("CaptureCardConfig_StartRecord Success\n", "green")
else
    WriteToLog("CaptureCardConfig_StartRecord failure\n", "red")
end
RSTD.Sleep(1000)

-- Start radar transmission
ar1.StartFrame()
WriteToLog("Streaming started — run updatedsar.py now\n", "green")
WriteToLog("Press Stop in mmWave Studio to end the session\n", "green")

-- Keep script alive
local elapsed = 0
while true do
    RSTD.Sleep(5000)
    elapsed = elapsed + 5
    WriteToLog("Streaming... " .. elapsed .. "s\n", "green")
end
-- Start radar transmission
ar1.StartFrame()
WriteToLog("Streaming started — run updatedsar.py now\n", "green")
WriteToLog("Press Stop in mmWave Studio to end the session\n", "green")

-- Keep script alive so radar keeps transmitting
local elapsed = 0
while true do
    RSTD.Sleep(5000)
    elapsed = elapsed + 5
    WriteToLog("Streaming... " .. elapsed .. "s\n", "green")
end