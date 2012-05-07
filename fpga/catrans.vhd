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

entity catrans is
    Port(
        CLK, R : in std_logic;
        I : in std_logic;
        N : in std_logic_vector(7 downto 0);
        O : out std_logic
    );
end catrans;

architecture Behavioral of catrans is
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
    
    component mux2 is
        generic(
            width : natural := 16
        );
        
        Port(
            Sel : in std_logic;
            I0, I1 : in std_logic_vector(width-1 downto 0);
            O : out std_logic_vector(width-1 downto 0)
        );
    end component;
    
    signal bRegD, bRegQ, sRegD, sRegQ : std_logic_vector(8 downto 0);
    
    signal cLvl0 : std_logic_vector(7 downto 0);
    signal cLvl1 : std_logic_vector(5 downto 0);
    
    signal tMask : std_logic_vector(8 downto 0);
    signal nSum : std_logic_vector(3 downto 0);
begin
    
    bRegD <= b"000001000";  -- hardcode Life config : B3
    sRegD <= b"000001100";  -- hardcode Life config : S23
    
    -- Life-like config B..S.. is held in registers
    -- TODO: make configurable from the outside
    -- TODO: store config up the hierarchy and let synthesis duplicate
    -- if performance benefits from duplication...
    bReg : reg
        generic map (width => 9)
        port map(
            CLK => CLK,
            E => '1',
            R => R,
            D => bRegD,
            Q => bRegQ
        );
    
    sReg : reg
        generic map (width => 9)
        port map(
            CLK => CLK,
            E => '1',
            R => R,
            D => sRegD,
            Q => sRegQ
        );
    
    
    -- select transfer mask based on input cell state
    cfg : mux2
        generic map (width => 9)
        port map(
            Sel => I,
            I0 => bRegQ,
            I1 => sRegQ,
            O => tMask
        );
    
    -- count live neighbours
    -- TODO: optimize?
    --nSum <= N(0) + N(1) + N(2) + N(3) + N(4) + N(5) + N(6) + N(7);
    
    cCount0 : for idx in 0 to 3 generate
        cLvl0(2*idx+0) <= N(2*idx+0) xor N(2*idx+1);
        cLvl0(2*idx+1) <= N(2*idx+0) and N(2*idx+1);
    end generate;
    
    cCount1 : for idx in 0 to 1 generate
        cLvl1(3*idx+0) <= cLvl0(4*idx+0) xor cLvl0(4*idx+2);
        cLvl1(3*idx+1) <= (cLvl0(4*idx+0) and cLvl0(4*idx+2)) or -- OR equiv XOR by def of inputs
                          (cLvl0(4*idx+1) xor cLvl0(4*idx+3));
        cLvl1(3*idx+2) <= cLvl0(4*idx+1) and cLvl0(4*idx+3);
    end generate;
    
    nSum(0) <= cLvl1(0) xor cLvl1(3);
    nSum(1) <= (cLvl1(1) xor cLvl1(4)) xor
               (cLvl1(0) and cLvl1(3));
    nSum(2) <= (cLvl1(2) xor cLvl1(5)) or
               (cLvl1(1) and cLvl1(4)) or
               ((cLvl1(1) xor cLvl1(4)) and (cLvl1(0) and cLvl1(3)));
    nSum(3) <= cLvl1(2) and cLvl1(5);
    
    -- use sum of neighbours states to select output value from trans mask
    O <= tMask(to_integer(unsigned(nSum)));
    
end Behavioral;
