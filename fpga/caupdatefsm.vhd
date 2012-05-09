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

entity caupdatefsm is
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
end caupdatefsm;

architecture Behavioral of caupdatefsm is
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
    
    type StateType is (
        SWait,   -- wait for coords : 00
        SFetchL, -- fetch x-1       : 01
        SFetchC, -- fetch x         : 11
        SFetchR  -- fetch x+1       : 10
    );
    
    -- state machine
    -- Manual state encoding for speed and ease of debugging
    signal sState, sNextState, sBaseState : std_logic_vector(1 downto 0);
    
    signal Ebase, Estate : std_logic;
begin
    cBaseState : reg
        generic map(width=> 2)
        port map(
            CLK=>CLK,
            E=>Ebase,
            R=>R,
            D=>sNextState,
            Q=>sBaseState
        );
    
    cState : reg
        generic map(width=> 2)
        port map(
            CLK=>CLK,
            E=>'1',
            R=>R,
            D=>sNextState,
            Q=>sState
        );
    
    -- state output
    process (sState, E, DRDYI, LASTY)
    begin
        -- default value to avoid latch inference and code duplication
        
        sNextState <= "00";
        
        case sState is
            when "00" =>  -- SWait
                if (E='1') then
                    -- determine cache validity (sBaseFetchState)
                    if (YHIT='1' and DISCARD='0' and (XHIT(0) xor XHIT(1))='1') then
                        if XHIT(0)='1' then
                            sNextState <= "10";  -- SFetchR
                        else
                            sNextState <= "11";  -- SFetchC
                        end if;
                    elsif (XBND='1') then
                        sNextState <= "11"; -- SFetchC;
                    else
                        sNextState <= "01"; -- SFetchL;
                    end if;
                end if;
                
            when "01" =>  -- SFetchL
                if (DRDYI='1') then
                    sNextState <= "11"; -- SFetchC;
                end if;
                
            when "11" => -- SFetchC
                if (DRDYI='1') then
                    sNextState <= "10"; -- sFetchR;
                end if;
                
            when "10" => -- SFetchR
                if (DRDYI='1') then
                    if (LASTY='1') then
                        -- macrocell finished
                        sNextState <= "00";
                    else
                        sNextState <= sBaseState;
                    end if;
                end if;
            when others =>
                sNextState <= "ZZ";
        end case;
    end process;
    
    Ebase <= E and not (sState(1) or sState(0));
    
    ELP <= Ebase;
    ELD <= DRDYI and (sState(1) or sState(0));
    CCC <= Ebase and XBND;
    CCR <= Ebase and YBND;
    DREQI <= sNextState(1) or sNextState(0);
    RDY <= not sNextState(1) and not sNextState(0);
    NXTX <= not sNextState(1) & not sNextState(0);
    NXTY <= (Ebase and YBND) or
            (sState(1) and not sState(0) and DRDYI);
end Behavioral;
