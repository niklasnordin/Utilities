//
//  functions.m
//  propertyCalculator
//
//  Created by Niklas Nordin on 2011-01-24.
//  Copyright 2011 nequam. All rights reserved.
//

#import "functions.h"


@implementation functions
 
// split the string between the spaces and return the first string
 -(NSString *) keyWord:(NSString *) line
 {
	 NSArray *words = [line componentsSeparatedByString:@" "];
	 return [words objectAtIndex:0];
 };
 
@end
