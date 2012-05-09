--------------------------------------------------------------------------------
-- Copyright 2012 Hugues Bruant <hugues.bruant@gmail.com>
-- All rights reserved.
--
-- This file is part of a school project and licensed under the terms of FreeBSD
-- license (2-clause BSD also refered to as Simplified BSD License)
--------------------------------------------------------------------------------

library ieee;
use ieee.std_logic_1164.all;

-- CA transition at macrocell level (restricted to Life-like for now)
-- IE : input enable
-- LC,RC : Left and right boundary condition respectively
-- I : macrocell state row
-- O : next macrocell state row
--
-- 2 clock cycle setup
-- 1 row/cycle in steady state [1 cycle latency when all neighbours are avail]
entity castream is
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
end castream;

architecture Behavioral of castream is
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
    
    component catrans is
        Port(
            CLK, R : in std_logic;
            I : in std_logic;
            N : in std_logic_vector(7 downto 0);
            O : out std_logic
        );
    end component;
    
    type regBarray is array (integer range <>) of std_logic_vector(width+1 downto 0);
    type regWarray is array (integer range <>) of std_logic_vector(width-1 downto 0);
    
    signal lvl : regBarray(3 downto 0);
    signal mix : regWarray(7 downto 0);
begin
    lvl(0) <= RC & I & LC;
    
    cMid : reg
        generic map(width => width+2)
        port map(
            CLK => CLK,
            E => IE,
            R => R,
            D => lvl(0),
            Q => lvl(1)
        );
    
    cBot : reg
        generic map(width => width+2)
        port map(
            CLK => CLK,
            E => IE,
            R => R,
            D => lvl(1),
            Q => lvl(2)
        );
    
    cTrans: for idx in 0 to width-1 generate
        -- neighbour mixing...
        mix(idx) <=
            lvl(0)(idx) & lvl(0)(idx+1) & lvl(0)(idx+2) &
            lvl(1)(idx) &                 lvl(1)(idx+2) &
            lvl(2)(idx) & lvl(2)(idx+1) & lvl(2)(idx+2);
        
        -- cell transition
        ctu : catrans
            port map(
                CLK => CLK,
                R => R,
                I => lvl(1)(idx+1),
                N => mix(idx),
                O => O(idx)
            );
    end generate;
end Behavioral;
