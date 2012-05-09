--------------------------------------------------------------------------------
-- Copyright 2012 Hugues Bruant <hugues.bruant@gmail.com>
-- All rights reserved.
--
-- This file is part of a school project and licensed under the terms of FreeBSD
-- license (2-clause BSD also refered to as Simplified BSD License)
--------------------------------------------------------------------------------

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

-- Macrocell-level CA update :
--  takes X,Y macrocell coords
--  fetch mem
--  update cell state
--  store new state
--
-- 
--
-- variable latency (depends on memory subsystem and cache-friendliness of
-- access pattern). assuming 1cc mem latency :
--  * full hit    : 18cc (Y unchanged, X inc by 1 from prev coords)
--  * partial hit : 36cc (Y unchanged, X inc by 2 from prev coords)
--  * full miss   : 54cc
--
-- steady state throughput : 32 cells/1cc w/ cache hit
--                           32 cells/3cc w/ cache miss

entity caupdate is
    generic(
        ybwidth : natural := 10;     -- 1024 rows
        xbwidth : natural := 5;      -- 1024 cols / 32 bits -> 32 macrocells
        dbwidth : natural := 5       -- 32 bits memory bus (macrocell width)
    );

    Port(
        CLK, E, R : in std_logic;
        
        -- do not use cache
        DISCARD : in std_logic;
        
        -- interface with work dispatcher
        X : in std_logic_vector(xbwidth-1 downto 0);
        Y : in std_logic_vector(ybwidth-1 downto 0);
        
        -- macrocell finished
        RDY : out std_logic;
        
        -- interface with memory subsystem
        ADDRI, ADDRO : out std_logic_vector((xbwidth + ybwidth)-1 downto 0);
        DREQI, DREQO : out std_logic;
        DRDYI, DRDYO : in std_logic;
        DATAI : in std_logic_vector((2**dbwidth)-1 downto 0);
        DATAO : out std_logic_vector((2**dbwidth)-1 downto 0)
    );
end caupdate;

architecture Behavioral of caupdate is
    constant streamwidth : natural := 8;
    constant logheight : natural := 4;     -- macrocell height
    constant height : natural := 2**logheight;
    constant datawidth : natural := 2**dbwidth;
    constant addrwidth : natural := xbwidth + ybwidth;
    constant numstream : natural := datawidth / streamwidth;
    
    component reg1 is
        Port(
            CLK, E, R : in std_logic;
            D : in std_logic;
            Q : out std_logic
        );
    end component;
    
    component reg is
        generic(
            width : natural := 16
        );
        
        Port(
            CLK, E, R : in std_logic;
            D : in std_logic_vector(width-1 downto 0);
            Q : out std_logic_vector(width-1 downto 0)
        );
    end component;
    
    component caupdatefsm is
        Port(
            CLK, E, R : in std_logic;
            
            -- do not use cache
            DISCARD : in std_logic;
            
            -- macrocell finished
            RDY : out std_logic;
            
            DREQI : out std_logic;
            DRDYI : in std_logic;
            
            --
            XBND, YBND : in std_logic;
            
            YHIT : in std_logic;
            XHIT : in std_logic_vector(1 downto 0);
            
            -- signals to datapath
            CCC, CCR : out std_logic; -- CacheClear Col/Row
            
            ELP : out std_logic;  -- pos load (start new macrocell)
            ELD : out std_logic;  -- data load
            
            NXTX : out std_logic_vector(1 downto 0);   -- L, C, R
            NXTY : out std_logic;               -- 0,+1 (start new row)
            LASTY : in std_logic
        );
    end component;
    
    component castream is
        generic(
            width : natural := 16
        );
        
        Port(
            CLK, R : in std_logic;
            
            IE : in std_logic;
            
            LC, RC : in std_logic;
            I : in std_logic_vector(width-1 downto 0);
            O : out std_logic_vector(width-1 downto 0)
        );
    end component;
    
    type regWarray is array (integer range <>) of std_logic_vector(datawidth-1 downto 0);
    
    -- FSM signals
    signal XBND, YBND : std_logic;
    signal YHIT : std_logic;
    signal XHIT : std_logic_vector(1 downto 0);
    signal CCC, CCR : std_logic; -- CacheClear Col/Row
    signal ELP, ELD : std_logic;  -- data load
    signal NXTX : std_logic_vector(1 downto 0);   -- L, C, R
    signal NXTY : std_logic;               -- 0,+1 (start new row)
    signal LASTY : std_logic;
    
    -- cache reset lines
    signal RcacheCur, RcacheNext : std_logic_vector(height+1 downto 0);
    
    -- register enable flags
    signal Enext, Efm : std_logic;
    signal Ecld : std_logic_vector(height+1 downto 0);
    
    -- helper for cache management
    signal diffY : std_logic_vector(ybwidth-1 downto 0);
    signal diffX : std_logic_vector(xbwidth-1 downto 0);
    
    -- row mask
    signal nfetchmask, rfetchmask : std_logic_vector(height+1 downto 0);
    signal rowidx, nrowidx : std_logic_vector(logheight downto 0);
    
    -- coords : cache, fetch position, ...
    signal cX, sfX, baseX : std_logic_vector(xbwidth-1 downto 0);
    signal cY, sfY, nfY, ssY, firstY, baseY, offY : std_logic_vector(ybwidth-1 downto 0);
    
    -- cached cell contents
    signal ccell, ncell : regWarray(height+1 downto 0);
    
    -- macrocell row to send to stream unit (w/ boundary)
    signal crow : std_logic_vector(datawidth+1 downto 0);
    
    signal sLastY, Rridx, sRDY, sRDYr, sreqo : std_logic;
begin
    -- cell location cache
    ccX : reg
        generic map(width => xbwidth)
        port map(
            CLK=>CLK,
            R=>R,
            E=>ELP,
            D=>X,
            Q=>cX
        );
    
    ccY : reg
        generic map(width => ybwidth)
        port map(
            CLK=>CLK,
            R=>R,
            E=>ELP,
            D=>Y,
            Q=>cY
        );
    
    -- register enable flag
    Enext <= NXTY;
    
    -- fetch coordinate
    cfY : reg
        generic map(width => ybwidth)
        port map(
            CLK=>CLK,
            R=>R,
            E=>Enext,
            D=>nfY,
            Q=>sfY
        );
    
    -- store coordinate : lagging behind fetch
    csY : reg
        generic map(width => ybwidth)
        port map(
            CLK=>CLK,
            R=>R,
            E=>Enext,
            D=>sfY,
            Q=>ssY
        );
    
    -- fetchmask
    nfetchmask <= (height downto 0 => '0') & '1' when R='1' else
                  rfetchmask(height downto 0) & rfetchmask(height+1);
    
    rdyEdge : reg1
        port map(
            CLK=>CLK,
            R=>'0',
            E=>'1',
            D=>sRDY,
            Q=>sRDYr
        );
    
    Efm <= NXTY or R or (sRDY and not sRDYr);
    
    cfm : reg
        generic map(width => height+2)
        port map(
            CLK => CLK,
            R => '0',
            E => Efm,
            D => nfetchmask,
            Q => rfetchmask
        );
    
    -- row counter
    Rridx <= R or sRDY;
    nrowidx <= std_logic_vector(unsigned(rowidx)+1);
    
    crowidx : reg
        generic map(width => 5)
        port map(
            CLK=>CLK,
            R=>Rridx,
            E=>NXTY,
            D=>nrowidx,
            Q=>rowidx
        );
    
    -- helper to compute cache validity
    diffY <= Y xor cY;
    diffX <= std_logic_vector(unsigned(x) - unsigned(cX));
    
    XHIT <= diffX(1 downto 0) when diffX(xbwidth-1 downto 2)=(xbwidth-3 downto 0=>'0') else "11";
    YHIT <= '1' when diffY=(ybwidth-1 downto 0=>'0') else '0';
    YBND <= '1' when Y=(ybwidth-1 downto 0=>'0') else '0';
    XBND <= '1' when X=(xbwidth-1 downto 0=>'0') else '0';
    
    sLastY <= '1' when rowidx="10001" else '0';
    LASTY <= sLastY;
    
    -- fetch address
    -- TODO: avoid overrun...
    baseX <= X when ELP='1' else cX;
    
    with NXTX select
        sfX <= std_logic_vector(unsigned(baseX) - 1) when "10",
               baseX                                 when "00",
               std_logic_vector(unsigned(baseX) + 1) when "01",
               (others => 'Z')                       when others;
    
    firstY <= Y when YBND='1' else std_logic_vector(unsigned(Y)-1);
    baseY <= firstY when ELP='1' else sfY;
    offY <= (offY'length-2 downto 0 => '0') & not ELP;
    nfY <= std_logic_vector(unsigned(baseY)+unsigned(offY));
    
    ADDRI(addrwidth-1 downto xbwidth) <= nfY when Enext='1' else baseY;
    ADDRI(xbwidth-1 downto 0) <= sfX;
    
    -- cell data cache
    
    RcacheCur(0) <= R or CCR or CCC;
    RcacheCur(height+1 downto 1) <= (others => R or CCC);
    
    cCurrentCell : for idx in 0 to height+1 generate
        Ecld(idx) <= ELD and rfetchmask(idx);
        
        ccc : reg
            generic map(width => datawidth)
            port map(
                CLK => CLK,
                R => RcacheCur(idx),
                E => Ecld(idx),
                D => ncell(idx),
                Q => ccell(idx)
            );
    end generate;
    
    RcacheNext(0) <= R or CCR;
    RcacheNext(height+1 downto 1) <= (others => R);
    
    cNextCell : for idx in 0 to height+1 generate
        cnc : reg
            generic map(width => datawidth)
            port map(
                CLK => CLK,
                R => RcacheNext(idx),
                E => Ecld(idx),
                D => DATAI,
                Q => ncell(idx)
            );
    end generate;
    
    -- FSM to control pipelined fetch, update and store of macrocell rows
    cFSM : caupdatefsm
        port map(
            CLK=>CLK,
            R=>R,
            E=>E,
            DISCARD=>DISCARD,
            RDY=>sRDY,
            DREQI=>DREQI,
            DRDYI=>DRDYI,
            XBND=>XBND,
            YBND=>YBND,
            YHIT=>YHIT,
            XHIT=>XHIT,
            CCC=>CCC,
            CCR=>CCR,
            ELP=>ELP,
            ELD=>ELD,
            NXTX=>NXTX,   -- L, C, R
            NXTY=>NXTY,   -- 0,+1 (start new row)
            LASTY=>LASTY
        );
    
    RDY <= sRDY;
    
    -- data mixing
    
    with rowidx select
        crow <=
            DATAI(0) & ncell( 0) & ccell( 0)(datawidth-1) when "00000",
            DATAI(0) & ncell( 1) & ccell( 1)(datawidth-1) when "00001",
            DATAI(0) & ncell( 2) & ccell( 2)(datawidth-1) when "00010",
            DATAI(0) & ncell( 3) & ccell( 3)(datawidth-1) when "00011",
            DATAI(0) & ncell( 4) & ccell( 4)(datawidth-1) when "00100",
            DATAI(0) & ncell( 5) & ccell( 5)(datawidth-1) when "00101",
            DATAI(0) & ncell( 6) & ccell( 6)(datawidth-1) when "00110",
            DATAI(0) & ncell( 7) & ccell( 7)(datawidth-1) when "00111",
            DATAI(0) & ncell( 8) & ccell( 8)(datawidth-1) when "01000",
            DATAI(0) & ncell( 9) & ccell( 9)(datawidth-1) when "01001",
            DATAI(0) & ncell(10) & ccell(10)(datawidth-1) when "01010",
            DATAI(0) & ncell(11) & ccell(11)(datawidth-1) when "01011",
            DATAI(0) & ncell(12) & ccell(12)(datawidth-1) when "01100",
            DATAI(0) & ncell(13) & ccell(13)(datawidth-1) when "01101",
            DATAI(0) & ncell(14) & ccell(14)(datawidth-1) when "01110",
            DATAI(0) & ncell(15) & ccell(15)(datawidth-1) when "01111",
            DATAI(0) & ncell(16) & ccell(16)(datawidth-1) when "10000",
            DATAI(0) & ncell(17) & ccell(17)(datawidth-1) when "10001",
            (others => 'Z') when others;
    
    -- cell transition stream
    cStream: for idx in 0 to numstream-1 generate
        csu : castream
            generic map(width=>streamwidth)
            port map(
                CLK=>CLK,
                R=>R,
                IE=>NXTY,
                LC=>crow(idx*streamwidth),
                RC=>crow((idx+1)*streamwidth+1),
                I=>crow((idx+1)*streamwidth downto idx*streamwidth+1),
                O=>DATAO((idx+1)*streamwidth-1 downto idx*streamwidth)
            );
    end generate;
    
    -- output results
    -- TODO: throttle input if output takes more than 1cc
    
    DREQO <= (NXTY or sRDY) and not rfetchmask(0) and not rfetchmask(1);
    
    ADDRO <= ssY & cX;
    
end Behavioral;
