//! Discovery only: explicit public certificates on conditional list pools.
use serde_json::{json,Value};
use std::collections::HashSet;
use std::io::{self,BufRead};
use sudoku_core::{CandidateState,GlobalColoringProof,GlobalColoringReplay};

fn peer(a:usize,b:usize)->bool { a!=b && (a/9==b/9 || a%9==b%9 || (a/27==b/27 && a%9/3==b%9/3)) }
fn component(pool:&[usize],target:usize)->Vec<usize> {
    let mut seen=vec![target]; let mut i=0;
    while i<seen.len() { for &c in pool { if !seen.contains(&c) && peer(seen[i],c) { seen.push(c); } } i+=1; }
    seen.sort();seen
}
fn main() {
 for line in io::stdin().lock().lines() {
    let case:Value=serde_json::from_str(&line.unwrap()).unwrap();
    let state:CandidateState=serde_json::from_value(case["state"].clone()).unwrap();
    let grid=state.to_grid().unwrap();
    let ds:Vec<u16>=serde_json::from_value(case["state"]["domains"].clone()).unwrap();
    let vals:Vec<u8>=serde_json::from_value(case["state"]["values"].clone()).unwrap();
    let sol=case["solution"].as_str().unwrap().as_bytes();
    let mut successes=vec![]; let mut visited=HashSet::new();let mut conclusions=HashSet::new(); let mut tested=0;
    let budget=case["pattern_budget"].as_u64().unwrap_or(12000) as usize;
    'palettes: for palette in (1u16..512).filter(|p|p.count_ones()==4 || p.count_ones()==3) {
      let pure:Vec<usize>=(0..81).filter(|&c|vals[c]==0 && ds[c]&!palette==0).collect();
      if pure.len()<8 {continue;}
      let mixed:Vec<usize>=(0..81).filter(|&c|vals[c]==0 && ds[c]&!palette!=0 && (ds[c]&palette).count_ones()>=2 && (ds[c]&!palette).count_ones()<=2).collect();
      for &target in &mixed {
        let extra=ds[target]&!palette;
        let mut endpoints=vec![];
        if extra.count_ones()==1 && sol[target]-b'0'==extra.trailing_zeros() as u8+1 { endpoints.push((extra.trailing_zeros() as u8+1,true)); }
        for d in 1..=9 {if ds[target]&(1<<(d-1))!=0 && sol[target]-b'0'!=d {endpoints.push((d,false));}}
        for (digit,value) in endpoints {
          if conclusions.contains(&(target,digit,value)) {continue;}
          for other in std::iter::once(None).chain(mixed.iter().filter(|&&c|c!=target).map(|&c|Some(c))) {
            let mut pool=pure.clone();pool.push(target); if let Some(c)=other {pool.push(c);}
            let cells=component(&pool,target);
            if cells.len()<8 || cells.len()>32 {continue;}
            if !visited.insert((cells.clone(),target,digit,value)) {continue;}
            tested+=1;if tested>budget {break 'palettes;}
            if let Some(proof)=GlobalColoringProof::from_cells(&grid,cells.clone(),target,digit,value) {
              assert!(proof.verify(&grid)); let replay=GlobalColoringReplay::capture(&grid,&proof).unwrap();assert!(replay.verify());
              successes.push(json!({"palette":palette,"graph_cells":cells,"conclusion":{"cell":target,"digit":digit,"value":value},"replay":replay,"extra_guardian":other}));
              conclusions.insert((target,digit,value));break;
            }
          }
        }
      }
    }
    println!("{}",json!({"name":case["name"],"patterns_tested":tested,"budget_exhausted":tested>budget,"successes":successes}));
 }
}
